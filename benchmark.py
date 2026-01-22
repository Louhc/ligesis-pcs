#!/usr/bin/env python3
"""
PCS Benchmark (Remote Server)

Supported schemes:
  - Single-thread: LigeSIS, DeepFold, Ligero
  - Distributed: dLigeSIS, dDeepFold, dDeepFoldBatch, etc.

Usage:
    # Interactive mode
    python3 benchmark.py

    # Single command
    python3 benchmark.py status
    python3 benchmark.py set-n 4                    # Set num_party=4
    python3 benchmark.py run -s ligesis -m 24       # Single-thread test
    python3 benchmark.py run -s dligesis -m 28      # Distributed test
"""

import argparse
import json
import math
import os
import re
import readline  # Command line history and editing support
import shlex
import subprocess
import sys
import tempfile
import threading
import time
from dataclasses import dataclass, asdict
from datetime import datetime
from pathlib import Path
from typing import Optional

# ============== Configuration ==============

WORKSPACE = Path(__file__).parent.resolve()
SERVERS_CONFIG = WORKSPACE / "ligesis-pcs" / "dTests" / "servers_16.json"
RESULTS_DIR = WORKSPACE / "bench_results"
ZONE = "us-central1-a"
REMOTE_DIR = "~/ligesis-pcs"

# Current number of parties
NUM_PARTY = 4

# Single-thread schemes
SINGLE_SCHEMES = {
    "ligesis": {"bench_name": "ligesis_bench", "display_name": "LigeSIS"},
    "deepfold": {"bench_name": "deepfold_bench", "display_name": "DeepFold"},
    "ligero": {"bench_name": "ligero_bench", "display_name": "Ligero"},
}

# Distributed schemes (example names)
DISTRIBUTED_SCHEMES = {
    "dligesis": {"example_name": "dLigesis", "display_name": "dLigeSIS"},
    "ddeepfold": {"example_name": "dDeepFold", "display_name": "dDeepFold"},
    "ddeepfoldbatch": {"example_name": "dDeepFoldBatch", "display_name": "dDeepFoldBatch"},
    "dmerkle": {"example_name": "dMerkle", "display_name": "dMerkle"},
    "dchunkedbatch": {"example_name": "dChunkedBatch", "display_name": "dChunkedBatch"},
    "dmultichunkedbatchbench": {"example_name": "dMultiChunkedBatchBench", "display_name": "dMultiChunkedBatchBench"},
    "dmultichunkedbatchprofile": {"example_name": "dMultiChunkedBatchProfile", "display_name": "dMultiChunkedBatchProfile"},
}

ALL_SCHEMES = {**SINGLE_SCHEMES, **DISTRIBUTED_SCHEMES}

DEFAULT_MUS = [24, 26, 28, 30]
DEFAULT_SINGLE_SCHEMES = ["ligesis", "deepfold", "ligero"]
DEFAULT_ITERATIONS = 1

# ============== Data Structures ==============

@dataclass
class BenchResult:
    scheme: str
    mu: int
    iteration: int
    timestamp: str
    success: bool
    num_parties: int = 1
    setup_time_ms: Optional[float] = None
    commit_time_ms: Optional[float] = None
    open_time_ms: Optional[float] = None
    verify_time_ms: Optional[float] = None
    prover_time_ms: Optional[float] = None
    total_time_ms: Optional[float] = None
    communication_bytes: Optional[int] = None
    raw_output: str = ""
    error: str = ""


# ============== Server Management ==============

def load_servers_config():
    with open(SERVERS_CONFIG) as f:
        return json.load(f)


def get_all_servers() -> list[dict]:
    config = load_servers_config()
    return config["servers"]


def get_active_servers() -> list[dict]:
    """Get server list for current NUM_PARTY"""
    all_servers = get_all_servers()
    return all_servers[:NUM_PARTY]


def get_server_name() -> str:
    """Get first server name (for single-thread tests)"""
    return get_active_servers()[0]["name"]


def get_user() -> str:
    config = load_servers_config()
    return config.get("user", "")


def get_zone() -> str:
    config = load_servers_config()
    return config.get("zone", ZONE)


def run_gcloud(args: list[str], timeout: int = 300) -> subprocess.CompletedProcess:
    cmd = ["gcloud"] + args
    return subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)


def gcloud_ssh(instance: str, command: str, timeout: int = 3600, retries: int = 3) -> dict:
    """Execute command via gcloud ssh with IAP tunneling support"""
    user = get_user()
    zone = get_zone()

    cmd = ["gcloud", "compute", "ssh"]
    if user:
        cmd.append(f"{user}@{instance}")
    else:
        cmd.append(instance)

    cmd.extend([
        "--zone", zone,
        "--tunnel-through-iap",  # Explicitly use IAP tunneling
        "--quiet",  # Reduce output noise
    ])

    # Add SSH options
    cmd.extend([
        "--ssh-flag=-o ServerAliveInterval=30",
        "--ssh-flag=-o ServerAliveCountMax=3",
        "--ssh-flag=-o ConnectTimeout=30",
        "--", "-T", command
    ])

    last_error = ""
    for attempt in range(retries):
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
            stderr = result.stderr

            # Filter out IAP tunneling warnings and other noise
            stderr_lines = [
                line for line in stderr.split('\n')
                if not line.startswith("External IP address was not found")
                and not line.startswith("WARNING:")
                and "Connection closed by UNKNOWN" not in line
                and "increase the performance of the tunnel" not in line
                and "installing NumPy" not in line
                and "cloud.google.com/iap/docs" not in line
                and line.strip()
            ]
            filtered_stderr = '\n'.join(stderr_lines)

            # Check for IAP connection error (needs retry)
            if result.returncode != 0 and "Connection closed by UNKNOWN" in stderr:
                last_error = f"IAP connection failed (attempt {attempt + 1}/{retries})"
                if attempt < retries - 1:
                    time.sleep(2 * (attempt + 1))  # Exponential backoff
                    continue

            return {
                "stdout": result.stdout,
                "stderr": filtered_stderr,
                "returncode": result.returncode
            }
        except subprocess.TimeoutExpired:
            last_error = "Timeout"
            if attempt < retries - 1:
                continue
            return {"stdout": "", "stderr": "Timeout", "returncode": -1}
        except Exception as e:
            last_error = str(e)
            if attempt < retries - 1:
                time.sleep(1)
                continue
            return {"stdout": "", "stderr": str(e), "returncode": -1}

    return {"stdout": "", "stderr": last_error, "returncode": -1}


def gcloud_scp(local_path: str, instance: str, remote_path: str, retries: int = 3) -> bool:
    user = get_user()
    zone = get_zone()

    cmd = [
        "gcloud", "compute", "scp",
        "--zone", zone,
        "--tunnel-through-iap",  # Explicitly use IAP
        "--quiet",
    ]
    remote_spec = f"{user}@{instance}:{remote_path}" if user else f"{instance}:{remote_path}"
    cmd.extend([local_path, remote_spec])

    for attempt in range(retries):
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
            if result.returncode == 0:
                return True
            # Retry on IAP connection error
            if "Connection closed by UNKNOWN" in result.stderr and attempt < retries - 1:
                time.sleep(2 * (attempt + 1))
                continue
            return False
        except subprocess.TimeoutExpired:
            if attempt < retries - 1:
                continue
            return False
        except Exception:
            if attempt < retries - 1:
                time.sleep(1)
                continue
            return False
    return False


def get_running_servers() -> list[str]:
    """Get all running node-* servers"""
    result = run_gcloud([
        "compute", "instances", "list",
        "--filter=name~'^node-' AND status=RUNNING",
        "--format=value(name)"
    ])
    if result.returncode != 0:
        return []
    return [s for s in result.stdout.strip().split('\n') if s]


def get_stopped_servers() -> list[str]:
    """Get all stopped node-* servers"""
    result = run_gcloud([
        "compute", "instances", "list",
        "--filter=name~'^node-' AND status=TERMINATED",
        "--format=value(name)"
    ])
    if result.returncode != 0:
        return []
    return [s for s in result.stdout.strip().split('\n') if s]


def cmd_status(_args=None):
    global NUM_PARTY
    print(f"\nCurrent config: num_party = {NUM_PARTY}")
    print(f"Active servers: node-1 to node-{NUM_PARTY}\n")

    result = run_gcloud(["compute", "instances", "list", "--filter=name~'^node-'"])
    print(result.stdout)
    return 0


def cmd_set_n(args):
    global NUM_PARTY
    n = args.n

    # Verify n is power of 2
    if n < 1 or (n & (n - 1)) != 0:
        print(f"Error: num_party must be a power of 2, got {n}")
        return 1

    all_servers = get_all_servers()
    if n > len(all_servers):
        print(f"Error: num_party ({n}) exceeds available servers ({len(all_servers)})")
        return 1

    NUM_PARTY = n
    print(f"Set num_party = {NUM_PARTY}")
    print(f"  Active servers: node-1 to node-{NUM_PARTY}")

    # Show current server status
    running = get_running_servers()
    needed = [f"node-{i}" for i in range(1, NUM_PARTY + 1)]
    extra = [s for s in running if s not in needed]

    if extra:
        print(f"\nWarning: Extra servers running: {', '.join(extra)}")
        print("  Use 'start' command to adjust server state")

    return 0


def cmd_start(_args=None):
    global NUM_PARTY

    needed = [f"node-{i}" for i in range(1, NUM_PARTY + 1)]
    running = get_running_servers()

    # Servers to start
    to_start = [s for s in needed if s not in running]
    # Servers to stop (extra ones)
    to_stop = [s for s in running if s not in needed]

    zone = get_zone()

    if to_stop:
        print(f"Stopping extra servers: {', '.join(to_stop)}")
        result = run_gcloud(["compute", "instances", "stop"] + to_stop + ["--zone", zone], timeout=180)
        if result.returncode != 0:
            print(f"Stop failed: {result.stderr}", file=sys.stderr)

    if to_start:
        print(f"Starting servers: {', '.join(to_start)}")
        result = run_gcloud(["compute", "instances", "start"] + to_start + ["--zone", zone], timeout=180)
        if result.returncode != 0:
            print(f"Start failed: {result.stderr}", file=sys.stderr)
            return 1

        print("Waiting for servers to be ready...")
        time.sleep(15)

    if not to_start and not to_stop:
        print(f"Servers ready: {', '.join(needed)} already running")
    else:
        print(f"Servers ready: {', '.join(needed)}")

    return 0


def cmd_stop(_args=None):
    running = get_running_servers()

    if not running:
        print("No running servers")
        return 0

    zone = get_zone()
    print(f"Stopping {len(running)} servers: {', '.join(running)}")
    run_gcloud(["compute", "instances", "stop"] + running + ["--zone", zone], timeout=180)
    print("Servers stopped")
    return 0


def cmd_sync(_args=None):
    """Sync code to all active servers (parallel)"""
    servers = get_active_servers()
    print(f"Syncing code to {len(servers)} servers...", end="", flush=True)

    # Create tarball
    tar_path = f"/tmp/ligesis_sync_{datetime.now().strftime('%H%M%S')}.tar.gz"
    exclude_args = ["--exclude=target", "--exclude=.git", "--exclude=bench_results", "--exclude=.claude"]

    tar_cmd = ["tar", "czf", tar_path] + exclude_args + ["-C", str(WORKSPACE), "."]
    result = subprocess.run(tar_cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f" failed")
        print(f"Failed to create tarball: {result.stderr}")
        return 1

    # Parallel sync
    sync_results = {}

    def do_sync(idx: int, server: dict):
        instance = server["name"]
        # Upload
        if not gcloud_scp(tar_path, instance, "~/ligesis_sync.tar.gz"):
            sync_results[idx] = (False, "upload failed")
            return
        # Extract
        result = gcloud_ssh(
            instance,
            f"mkdir -p {REMOTE_DIR} && cd {REMOTE_DIR} && rm -rf * && "
            f"tar xzf ~/ligesis_sync.tar.gz && rm ~/ligesis_sync.tar.gz",
            timeout=120
        )
        if result["returncode"] != 0:
            sync_results[idx] = (False, f"extract failed: {result['stderr']}")
            return
        sync_results[idx] = (True, "ok")

    threads = []
    for i, server in enumerate(servers):
        t = threading.Thread(target=do_sync, args=(i, server))
        threads.append(t)

    for t in threads:
        t.start()
    for t in threads:
        t.join()

    Path(tar_path).unlink(missing_ok=True)

    # Check results
    failed = [(i, sync_results[i][1]) for i in range(len(servers)) if not sync_results.get(i, (False, "unknown"))[0]]
    if failed:
        print(" failed")
        for i, err in failed:
            print(f"  x {servers[i]['name']}: {err}")
        return 1

    print(" done")
    return 0


# ============== Parsers ==============

def parse_duration(s: str) -> Optional[float]:
    s = s.strip()
    patterns = [
        (r'([\d.]+)\s*s$', lambda x: float(x) * 1000),
        (r'([\d.]+)\s*ms$', lambda x: float(x)),
        (r'([\d.]+)\s*us$', lambda x: float(x) / 1000),
        (r'([\d.]+)\s*µs$', lambda x: float(x) / 1000),
        (r'([\d.]+)\s*ns$', lambda x: float(x) / 1_000_000),
    ]
    for pattern, converter in patterns:
        m = re.search(pattern, s, re.IGNORECASE)
        if m:
            return converter(m.group(1))
    return None


def parse_benchmark_output(output: str) -> dict:
    result = {}
    # Use [^,\n]+ to match until comma or newline (handles comma-separated format)
    # Priority: avg values first (for multi-iteration runs), then single values
    patterns = {
        "setup": [r'Setup \(avg\)[:\s]+([^,\n]+)', r'Setup[:\s]+([^,\n]+)'],
        "commit": [r'Commit \(avg\)[:\s]+([^,\n]+)', r'Commit[:\s]+([^,\n]+)'],
        "open": [r'Open \(avg\)[:\s]+([^,\n]+)', r'Open[:\s]+([^,\n]+)'],
        "verify": [r'Verify \(avg\)[:\s]+([^,\n]+)', r'Verify[:\s]+([^,\n]+)'],
        "total": [r'Total \(avg\)[:\s]+([^,\n]+)', r'Total[:\s]+([^,\n]+)'],
    }

    for key, pats in patterns.items():
        for pat in pats:
            m = re.search(pat, output, re.IGNORECASE)
            if m:
                duration = parse_duration(m.group(1))
                if duration is not None:
                    result[key] = duration
                    break

    # Parse communication stats (try BYTES first, then MB)
    comm_match = re.search(r'COMM_TOTAL_BYTES:\s*(\d+)', output)
    if comm_match:
        result["communication_bytes"] = int(comm_match.group(1))
    else:
        comm_mb_match = re.search(r'COMM_TOTAL_MB:\s*([\d.]+)', output)
        if comm_mb_match:
            result["communication_bytes"] = int(float(comm_mb_match.group(1)) * 1024 * 1024)

    return result


# ============== Single-thread Benchmark ==============

def check_bench_exists(scheme: str) -> bool:
    if scheme not in SINGLE_SCHEMES:
        return False
    bench_name = SINGLE_SCHEMES[scheme]["bench_name"]
    bench_path = WORKSPACE / "ligesis-pcs" / "benches" / f"{bench_name}.rs"
    return bench_path.exists()


def run_single_thread_benchmark(scheme: str, mu: int, iterations: int = 1) -> BenchResult:
    config = SINGLE_SCHEMES.get(scheme)
    if not config:
        return BenchResult(
            scheme=scheme, mu=mu, iteration=1,
            timestamp=datetime.now().isoformat(),
            success=False, error=f"Unknown scheme: {scheme}"
        )

    bench_name = config["bench_name"]
    instance = get_server_name()
    cmd = (
        f"cd {REMOTE_DIR} && source ~/.cargo/env && "
        f"cargo bench --package ligesis-pcs --bench {bench_name} "
        f"--features print-trace -- --mu {mu} --iterations {iterations} 2>&1"
    )

    print(f"  Running...", end="", flush=True)
    result = gcloud_ssh(instance, cmd, timeout=3600)
    output = result["stdout"] + result["stderr"]

    if result["returncode"] != 0:
        print(" failed")
        return BenchResult(
            scheme=scheme, mu=mu, iteration=iterations,
            timestamp=datetime.now().isoformat(),
            success=False, error=output[-500:], raw_output=output
        )

    parsed = parse_benchmark_output(output)
    if not parsed:
        print(" parse failed")
        return BenchResult(
            scheme=scheme, mu=mu, iteration=iterations,
            timestamp=datetime.now().isoformat(),
            success=False, error="Failed to parse output", raw_output=output
        )

    commit_ms = parsed.get("commit")
    open_ms = parsed.get("open")
    prover_ms = (commit_ms or 0) + (open_ms or 0) if commit_ms or open_ms else None

    print(" done")
    return BenchResult(
        scheme=scheme, mu=mu, iteration=iterations,
        timestamp=datetime.now().isoformat(),
        success=True,
        setup_time_ms=parsed.get("setup"),
        commit_time_ms=commit_ms,
        open_time_ms=open_ms,
        verify_time_ms=parsed.get("verify"),
        prover_time_ms=prover_ms,
        total_time_ms=parsed.get("total"),
        raw_output=output,
    )


# ============== Distributed Benchmark ==============

def compute_optimal_base_mu(mu: int, num_parties: int) -> int:
    log_parties = int(math.log2(num_parties))
    local_num_vars = mu - log_parties
    OPTIMAL_BASE_MU = 14
    return min(OPTIMAL_BASE_MU, local_num_vars)


def generate_network_config(hosts: list[str], base_port: int = 18000) -> str:
    """Generate network config (all nodes use same port)"""
    return "\n".join(f"{host}:{base_port}" for host in hosts)


def _run_gcloud_ssh_worker(instance: str, command: str, timeout: int, result_dict: dict, index: int):
    """Worker function for threaded gcloud ssh execution"""
    result_dict[index] = gcloud_ssh(instance, command, timeout=timeout)


def check_servers_running() -> tuple[bool, list[str], list[str]]:
    """Check if required servers are running, returns (all_running, running_list, not_running_list)"""
    needed = [f"node-{i}" for i in range(1, NUM_PARTY + 1)]
    running = get_running_servers()
    not_running = [s for s in needed if s not in running]
    return len(not_running) == 0, [s for s in needed if s in running], not_running


def run_distributed_benchmark(
    scheme: str,
    mu: int,
    iterations: int = 1,
    trace: bool = True,
    build: bool = False,
    base_mu: Optional[int] = None,
) -> BenchResult:
    global NUM_PARTY

    # Check server status
    all_running, running, not_running = check_servers_running()
    if not all_running:
        return BenchResult(
            scheme=scheme, mu=mu, iteration=iterations,
            timestamp=datetime.now().isoformat(),
            success=False, num_parties=NUM_PARTY,
            error=f"Servers not running: {', '.join(not_running)}. Run 'start' first"
        )

    config = DISTRIBUTED_SCHEMES.get(scheme)
    if not config:
        return BenchResult(
            scheme=scheme, mu=mu, iteration=iterations,
            timestamp=datetime.now().isoformat(),
            success=False, num_parties=NUM_PARTY,
            error=f"Unknown distributed scheme: {scheme}"
        )

    example_name = config["example_name"]
    servers = get_active_servers()
    num_parties = len(servers)

    # Compute base_mu
    actual_base_mu = base_mu if base_mu is not None else compute_optimal_base_mu(mu, num_parties)

    # Generate network config
    server_config = load_servers_config()
    hosts = [s["host"] for s in servers]
    network_port = server_config.get("network_port", 18000)
    config_content = generate_network_config(hosts, network_port)

    print(f"  Nodes: {num_parties}, base_mu: {actual_base_mu}")

    # Build if needed
    if build:
        print(f"  Building on remote servers...", end="", flush=True)
        build_cmd = (
            f"source ~/.cargo/env && cd {REMOTE_DIR}/ligesis-pcs && "
            f"RUSTFLAGS='-Awarnings' cargo build --example {example_name} --release"
        )
        if trace:
            build_cmd += " --features print-trace"
        build_cmd += " 2>&1"

        # Parallel build
        build_results = {}
        build_threads = []
        for i, server in enumerate(servers):
            t = threading.Thread(
                target=_run_gcloud_ssh_worker,
                args=(server["name"], build_cmd, 600, build_results, i)
            )
            build_threads.append(t)

        for t in build_threads:
            t.start()
        for t in build_threads:
            t.join()

        # Check build results
        build_failed = [i for i in range(num_parties) if build_results.get(i, {}).get("returncode", -1) != 0]
        if build_failed:
            print(" failed")
            return BenchResult(
                scheme=scheme, mu=mu, iteration=iterations,
                timestamp=datetime.now().isoformat(),
                success=False, num_parties=num_parties,
                error=f"Build failed on {', '.join(servers[i]['name'] for i in build_failed)}",
                raw_output=build_results.get(build_failed[0], {}).get("stdout", "")
            )
        print(" done")

    # Deploy network config (parallel)
    print(f"  Deploying network config...", end="", flush=True)
    config_cmd = f"cat > /tmp/ligesis_network.conf << 'EOF'\n{config_content}\nEOF"

    config_results = {}
    config_threads = []
    for i, server in enumerate(servers):
        t = threading.Thread(
            target=_run_gcloud_ssh_worker,
            args=(server["name"], config_cmd, 60, config_results, i)
        )
        config_threads.append(t)

    for t in config_threads:
        t.start()
    for t in config_threads:
        t.join()

    # Check config results
    config_failed = [i for i in range(num_parties) if config_results.get(i, {}).get("returncode", -1) != 0]
    if config_failed:
        print(" failed")
        return BenchResult(
            scheme=scheme, mu=mu, iteration=iterations,
            timestamp=datetime.now().isoformat(),
            success=False, num_parties=num_parties,
            error=f"Config deploy failed on {', '.join(servers[i]['name'] for i in config_failed)}"
        )
    print(" done")

    # Run test - needs concurrent execution (distributed protocol requires simultaneous run)
    print(f"  Running...", end="", flush=True)

    binary_path = f"{REMOTE_DIR}/target/release/examples/{example_name}"
    run_results = {}
    threads = []

    # Start all threads using separate worker function to avoid closure issues
    for i, server in enumerate(servers):
        run_cmd = (
            f"RUST_BACKTRACE=1 {binary_path} {i} /tmp/ligesis_network.conf "
            f"--mu {mu} --base-mu {actual_base_mu} --iterations {iterations} 2>&1"
        )
        t = threading.Thread(
            target=_run_gcloud_ssh_worker,
            args=(server["name"], run_cmd, 1800, run_results, i)
        )
        threads.append(t)

    # Start all threads
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    # Check results
    failed = [i for i in range(num_parties) if run_results.get(i, {}).get("returncode", -1) != 0]
    if failed:
        print(" failed")
        errors = []
        for i in failed:
            # Check both stdout and stderr for error messages
            stdout = run_results.get(i, {}).get("stdout", "")
            stderr = run_results.get(i, {}).get("stderr", "")
            # Look for actual errors in stdout (since we use 2>&1)
            err_lines = [l for l in stdout.split('\n') if 'error' in l.lower() or 'No such file' in l or 'not found' in l.lower()]
            if err_lines:
                err = err_lines[0]
            elif stderr:
                err = stderr
            else:
                err = stdout[-200:] if stdout else "Unknown error"
            errors.append(f"Party {i}: {err[:200]}")
        return BenchResult(
            scheme=scheme, mu=mu, iteration=iterations,
            timestamp=datetime.now().isoformat(),
            success=False, num_parties=num_parties,
            error="\n".join(errors),
            raw_output=run_results.get(0, {}).get("stdout", "")
        )

    # Parse output (from party 0)
    output = run_results[0]["stdout"]
    parsed = parse_benchmark_output(output)

    print(" done")

    commit_ms = parsed.get("commit")
    open_ms = parsed.get("open")
    prover_ms = (commit_ms or 0) + (open_ms or 0) if commit_ms or open_ms else None

    return BenchResult(
        scheme=scheme, mu=mu, iteration=iterations,
        timestamp=datetime.now().isoformat(),
        success=True, num_parties=num_parties,
        setup_time_ms=parsed.get("setup"),
        commit_time_ms=commit_ms,
        open_time_ms=open_ms,
        verify_time_ms=parsed.get("verify"),
        prover_time_ms=prover_ms,
        total_time_ms=parsed.get("total"),
        communication_bytes=parsed.get("communication_bytes"),
        raw_output=output,
    )


# ============== Result Management ==============

def save_result(result: BenchResult, batch_id: Optional[str] = None):
    if result.num_parties > 1:
        result_dir = RESULTS_DIR / f"distributed_n{result.num_parties}"
    else:
        result_dir = RESULTS_DIR / "single_thread"

    result_dir.mkdir(parents=True, exist_ok=True)

    if batch_id:
        result_file = result_dir / f"batch_{batch_id}.jsonl"
    else:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        result_file = result_dir / f"{result.scheme}_mu{result.mu}_{timestamp}.json"

    with open(result_file, 'a') as f:
        f.write(json.dumps(asdict(result)) + "\n")
    return result_file


def load_all_results(distributed: bool = False, num_parties: Optional[int] = None) -> list[BenchResult]:
    """Load all results from results directory"""
    results = []

    if distributed and num_parties:
        search_dirs = [RESULTS_DIR / f"distributed_n{num_parties}"]
    elif distributed:
        # Load all distributed results
        search_dirs = list(RESULTS_DIR.glob("distributed_n*"))
    else:
        search_dirs = [RESULTS_DIR / "single_thread"]

    for result_dir in search_dirs:
        if not result_dir.exists():
            continue
        for f in result_dir.glob("*.json*"):
            with open(f) as fp:
                for line in fp:
                    line = line.strip()
                    if line:
                        try:
                            data = json.loads(line)
                            results.append(BenchResult(**data))
                        except (json.JSONDecodeError, TypeError):
                            pass
    return results


def print_summary_table(results: list[BenchResult]):
    """Print summary table of results"""
    by_scheme = {}
    for r in results:
        by_scheme.setdefault(r.scheme, []).append(r)

    all_mus = sorted(set(r.mu for r in results))
    is_distributed = any(r.num_parties > 1 for r in results)

    if is_distributed:
        header = f"{'Scheme':<12} {'mu':>4} {'n':>3} | {'Commit':>10} {'Open':>10} {'Prover':>10} | {'Verify':>10} {'Comm':>10}"
    else:
        header = f"{'Scheme':<12} {'mu':>4} | {'Commit':>10} {'Open':>10} {'Prover':>10} | {'Verify':>10}"
    print(header)
    print("-" * len(header))

    for scheme in sorted(by_scheme.keys()):
        scheme_results = {r.mu: r for r in by_scheme[scheme]}
        display_name = ALL_SCHEMES.get(scheme, {}).get("display_name", scheme)
        for mu in all_mus:
            r = scheme_results.get(mu)
            if r:
                if is_distributed:
                    print(f"{display_name:<12} {mu:>4} {r.num_parties:>3} | "
                          f"{format_time(r.commit_time_ms):>10} "
                          f"{format_time(r.open_time_ms):>10} "
                          f"{format_time(r.prover_time_ms):>10} | "
                          f"{format_time(r.verify_time_ms):>10} "
                          f"{format_bytes(r.communication_bytes):>10}")
                else:
                    print(f"{display_name:<12} {mu:>4} | "
                          f"{format_time(r.commit_time_ms):>10} "
                          f"{format_time(r.open_time_ms):>10} "
                          f"{format_time(r.prover_time_ms):>10} | "
                          f"{format_time(r.verify_time_ms):>10}")


def export_csv(results: list[BenchResult], filename: str):
    """Export results to CSV file"""
    import csv
    is_distributed = any(r.num_parties > 1 for r in results)

    with open(filename, 'w', newline='') as f:
        writer = csv.writer(f)
        if is_distributed:
            writer.writerow(["Scheme", "mu", "Parties", "Commit (ms)", "Open (ms)", "Prover (ms)",
                           "Verify (ms)", "Total (ms)", "Comm (bytes)", "Timestamp"])
        else:
            writer.writerow(["Scheme", "mu", "Commit (ms)", "Open (ms)", "Prover (ms)",
                           "Verify (ms)", "Total (ms)", "Timestamp"])

        for r in results:
            display_name = ALL_SCHEMES.get(r.scheme, {}).get("display_name", r.scheme)
            row = [
                display_name,
                r.mu,
            ]
            if is_distributed:
                row.append(r.num_parties)
            row.extend([
                f"{r.commit_time_ms:.2f}" if r.commit_time_ms else "",
                f"{r.open_time_ms:.2f}" if r.open_time_ms else "",
                f"{r.prover_time_ms:.2f}" if r.prover_time_ms else "",
                f"{r.verify_time_ms:.2f}" if r.verify_time_ms else "",
                f"{r.total_time_ms:.2f}" if r.total_time_ms else "",
            ])
            if is_distributed:
                row.append(r.communication_bytes if r.communication_bytes else "")
            row.append(r.timestamp)
            writer.writerow(row)


def format_time(ms: Optional[float]) -> str:
    if ms is None:
        return "-"
    if ms >= 1000:
        return f"{ms/1000:.2f}s"
    return f"{ms:.1f}ms"


def format_bytes(b: Optional[int]) -> str:
    if b is None:
        return "-"
    if b >= 1024 * 1024:
        return f"{b / (1024 * 1024):.2f}MB"
    if b >= 1024:
        return f"{b / 1024:.2f}KB"
    return f"{b}B"


# ============== Command Handlers ==============

def cmd_run(args):
    scheme = args.scheme.lower()
    mu = args.mu
    iterations = args.iterations
    build = getattr(args, 'build', False)

    if scheme not in ALL_SCHEMES:
        print(f"Unknown scheme: {scheme}")
        print(f"Available schemes: {', '.join(ALL_SCHEMES.keys())}")
        return 1

    is_distributed = scheme in DISTRIBUTED_SCHEMES
    display_name = ALL_SCHEMES[scheme]['display_name']

    print(f"\n{'='*60}")
    if is_distributed:
        print(f"Running: {display_name}, mu={mu}, nodes={NUM_PARTY}, iter={iterations}")
    else:
        print(f"Running: {display_name}, mu={mu}")
    print(f"{'='*60}\n")

    if is_distributed:
        result = run_distributed_benchmark(scheme, mu, iterations, build=build)
    else:
        result = run_single_thread_benchmark(scheme, mu, iterations)

    if result.success:
        print(f"\nSuccess")
        print(f"  Setup:   {format_time(result.setup_time_ms)}")
        print(f"  Commit:  {format_time(result.commit_time_ms)}")
        print(f"  Open:    {format_time(result.open_time_ms)}")
        print(f"  Verify:  {format_time(result.verify_time_ms)}")
        print(f"  Prover:  {format_time(result.prover_time_ms)}")
        if result.communication_bytes:
            print(f"  Comm:    {format_bytes(result.communication_bytes)}")
    else:
        print(f"\nFailed: {result.error[:300]}")

    result_file = save_result(result)
    print(f"\nResult saved to: {result_file}")
    return 0 if result.success else 1


def cmd_batch(args):
    schemes = [s.strip().lower() for s in args.schemes.split(',')] if args.schemes else DEFAULT_SINGLE_SCHEMES
    mus = [int(m.strip()) for m in args.mus.split(',')] if args.mus else DEFAULT_MUS
    iterations = args.iterations

    available_schemes = []
    for scheme in schemes:
        if scheme not in ALL_SCHEMES:
            print(f"Warning: Unknown scheme: {scheme}, skipping")
            continue
        if scheme in SINGLE_SCHEMES and not check_bench_exists(scheme):
            print(f"Warning: {scheme} benchmark not found, skipping")
            continue
        available_schemes.append(scheme)

    if not available_schemes:
        print("No available benchmarks")
        return 1

    tests = [(scheme, mu) for scheme in available_schemes for mu in mus]
    total = len(tests)

    print(f"\n{'='*60}")
    print(f"PCS Benchmark Batch")
    print(f"Schemes: {', '.join(available_schemes)}")
    print(f"mu: {', '.join(map(str, mus))}")
    print(f"Iterations: {iterations}")
    print(f"Total tests: {total}")
    print(f"{'='*60}\n")

    batch_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    results = []

    for i, (scheme, mu) in enumerate(tests, 1):
        display_name = ALL_SCHEMES[scheme]['display_name']
        is_distributed = scheme in DISTRIBUTED_SCHEMES

        if is_distributed:
            print(f"[{i}/{total}] {display_name} mu={mu} n={NUM_PARTY}")
            result = run_distributed_benchmark(scheme, mu, iterations)
        else:
            print(f"[{i}/{total}] {display_name} mu={mu}")
            result = run_single_thread_benchmark(scheme, mu, iterations)

        results.append(result)
        save_result(result, batch_id)

        if result.success:
            extra = f", comm={format_bytes(result.communication_bytes)}" if result.communication_bytes else ""
            print(f"        prover={format_time(result.prover_time_ms)}, verify={format_time(result.verify_time_ms)}{extra}\n")
        else:
            print(f"        Error: {result.error[:80]}\n")

    print(f"\n{'='*60}")
    print("Batch complete!")
    print(f"{'='*60}\n")

    # Print result file locations
    if results:
        saved_dirs = set()
        for r in results:
            if r.num_parties > 1:
                saved_dirs.add(RESULTS_DIR / f"distributed_n{r.num_parties}")
            else:
                saved_dirs.add(RESULTS_DIR / "single_thread")
        for d in saved_dirs:
            result_file = d / f"batch_{batch_id}.jsonl"
            print(f"Results saved to: {result_file}")

    return 0


def cmd_report(args):
    """Show/export benchmark results"""
    distributed = getattr(args, 'distributed', False)
    num_parties = getattr(args, 'n', None)

    results = load_all_results(distributed=distributed, num_parties=num_parties)
    if not results:
        print("No results found")
        return 1

    successful = [r for r in results if r.success]
    if not successful:
        print("No successful test results")
        return 1

    # Get latest result for each (scheme, mu) pair
    latest = {}
    for r in sorted(successful, key=lambda x: x.timestamp):
        latest[(r.scheme, r.mu)] = r

    results_list = sorted(latest.values(), key=lambda x: (x.scheme, x.mu))
    print_summary_table(results_list)

    csv_file = getattr(args, 'csv', None)
    if csv_file:
        export_csv(results_list, csv_file)
        print(f"\nCSV exported to: {csv_file}")
    return 0


def cmd_list(_args=None):
    print("\nAvailable Benchmarks:")
    print("-" * 50)
    print("\nSingle-thread:")
    for scheme, config in SINGLE_SCHEMES.items():
        exists = "+" if check_bench_exists(scheme) else "-"
        print(f"  {exists} {config['display_name']:<20} ({scheme})")

    print("\nDistributed:")
    for scheme, config in DISTRIBUTED_SCHEMES.items():
        print(f"  + {config['display_name']:<20} ({scheme})")
    print()
    return 0


def cmd_help(_args=None):
    global NUM_PARTY
    print(f"""
Current config: num_party = {NUM_PARTY}

Commands:
  status              Show server status
  set-n <n>           Set num_party (must be power of 2)
  start               Start required servers (node-1 to node-{NUM_PARTY})
  stop                Stop all servers
  sync                Sync code to servers
  list                List available benchmarks

  run -s <scheme> -m <mu> [-i <iterations>] [--build]
                      Run benchmark
                      Example: run -s ligesis -m 24      # single-thread
                      Example: run -s dligesis -m 28     # distributed

  batch [-s <schemes>] [-m <mus>] [-i <iterations>]
                      Run batch benchmarks
                      Example: batch -s ligesis,deepfold -m 24,26,28

  report [--csv <file>] [-d] [-n <parties>]
                      Show/export results
                      -d, --distributed   Show distributed results
                      -n <parties>        Filter by num_parties
                      --csv <file>        Export to CSV file

  help                Show this help
  exit, quit, q       Exit
""")
    return 0


# ============== Interactive Mode ==============

def create_parser():
    parser = argparse.ArgumentParser(prog="", add_help=False)
    subparsers = parser.add_subparsers(dest="command")

    subparsers.add_parser("status")
    subparsers.add_parser("start")
    subparsers.add_parser("stop")
    subparsers.add_parser("sync")
    subparsers.add_parser("list")
    subparsers.add_parser("help")

    p = subparsers.add_parser("set-n")
    p.add_argument("n", type=int, help="Number of parties")

    p = subparsers.add_parser("run")
    p.add_argument("--scheme", "-s", type=str, required=True)
    p.add_argument("--mu", "-m", type=int, default=24)
    p.add_argument("--iterations", "-i", type=int, default=1)
    p.add_argument("--build", "-b", action="store_true", help="Build on remote before running")

    p = subparsers.add_parser("batch")
    p.add_argument("--schemes", "-s", type=str, default=None)
    p.add_argument("--mus", "-m", type=str, default=None)
    p.add_argument("--iterations", "-i", type=int, default=DEFAULT_ITERATIONS)

    p = subparsers.add_parser("report")
    p.add_argument("--csv", type=str, default=None, help="Export to CSV file")
    p.add_argument("--distributed", "-d", action="store_true", help="Show distributed results")
    p.add_argument("-n", type=int, default=None, help="Filter by num_parties")

    return parser


def interactive_mode():
    print("=" * 60)
    print("PCS Benchmark - Interactive Mode")
    print(f"Current num_party = {NUM_PARTY}")
    print("Type 'help' for commands, 'exit' to quit")
    print("=" * 60)

    parser = create_parser()
    cmd_map = {
        "status": cmd_status,
        "set-n": cmd_set_n,
        "start": cmd_start,
        "stop": cmd_stop,
        "sync": cmd_sync,
        "list": cmd_list,
        "help": cmd_help,
        "run": cmd_run,
        "batch": cmd_batch,
        "report": cmd_report,
    }

    while True:
        try:
            line = input("\n> ").strip()
        except (EOFError, KeyboardInterrupt):
            print("\nExit")
            break

        if not line:
            continue

        if line.lower() in ("exit", "quit", "q"):
            print("Exit")
            break

        try:
            argv = shlex.split(line)
            args = parser.parse_args(argv)

            if args.command in cmd_map:
                cmd_map[args.command](args)
            else:
                print(f"Unknown command: {line}")
                cmd_help()

        except SystemExit:
            pass
        except Exception as e:
            print(f"Error: {e}")


# ============== Main ==============

def main():
    global NUM_PARTY

    if len(sys.argv) == 1:
        interactive_mode()
        return 0

    parser = argparse.ArgumentParser(
        description="PCS Benchmark (Remote Server)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s                                         # Interactive mode
  %(prog)s status                                  # Show server status
  %(prog)s set-n 4                                 # Set num_party=4
  %(prog)s start                                   # Start servers
  %(prog)s sync                                    # Sync code
  %(prog)s run -s ligesis -m 24                    # Single-thread test
  %(prog)s run -s dligesis -m 28 --build           # Distributed test
        """
    )

    subparsers = parser.add_subparsers(dest="command", help="Commands")

    subparsers.add_parser("status", help="Show server status").set_defaults(func=cmd_status)
    subparsers.add_parser("start", help="Start servers").set_defaults(func=cmd_start)
    subparsers.add_parser("stop", help="Stop servers").set_defaults(func=cmd_stop)
    subparsers.add_parser("sync", help="Sync code to servers").set_defaults(func=cmd_sync)
    subparsers.add_parser("list", help="List available benchmarks").set_defaults(func=cmd_list)

    p = subparsers.add_parser("set-n", help="Set num_party")
    p.add_argument("n", type=int, help="Number of parties (must be power of 2)")
    p.set_defaults(func=cmd_set_n)

    p = subparsers.add_parser("run", help="Run benchmark")
    p.add_argument("--scheme", "-s", type=str, required=True)
    p.add_argument("--mu", "-m", type=int, default=24)
    p.add_argument("--iterations", "-i", type=int, default=1)
    p.add_argument("--build", "-b", action="store_true", help="Build on remote before running")
    p.set_defaults(func=cmd_run)

    p = subparsers.add_parser("batch", help="Run batch benchmarks")
    p.add_argument("--schemes", "-s", type=str, default=None,
                   help=f"Comma-separated schemes (default: {','.join(DEFAULT_SINGLE_SCHEMES)})")
    p.add_argument("--mus", "-m", type=str, default=None,
                   help=f"Comma-separated mus (default: {','.join(map(str, DEFAULT_MUS))})")
    p.add_argument("--iterations", "-i", type=int, default=DEFAULT_ITERATIONS)
    p.set_defaults(func=cmd_batch)

    p = subparsers.add_parser("report", help="Show/export results")
    p.add_argument("--csv", type=str, default=None, help="Export to CSV file")
    p.add_argument("--distributed", "-d", action="store_true", help="Show distributed results")
    p.add_argument("-n", type=int, default=None, help="Filter by num_parties")
    p.set_defaults(func=cmd_report)

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        return 1

    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
