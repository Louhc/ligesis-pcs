#!/usr/bin/env python3
"""
Single-Thread PCS Benchmark (Remote Server)

支持的方案: LigeSIS, DeepFold, Ligero

使用方法:
    # 交互式模式
    python3 single_thread_benchmark.py

    # 单次命令
    python3 single_thread_benchmark.py status
    python3 single_thread_benchmark.py run --scheme ligesis --mu 24
    python3 single_thread_benchmark.py batch --schemes ligesis,deepfold --mus 24,26,28 -i 3
"""

import argparse
import json
import re
import readline  # 支持命令行历史和编辑
import shlex
import subprocess
import sys
import time
from dataclasses import dataclass, asdict
from datetime import datetime
from pathlib import Path
from typing import Optional

# ============== 配置 ==============

WORKSPACE = Path(__file__).parent.resolve()
SERVERS_CONFIG = WORKSPACE / "ligesis-pcs" / "dTests" / "servers_16.json"
RESULTS_DIR = WORKSPACE / "bench_results" / "single_thread"
ZONE = "us-central1-a"

REMOTE_DIR = "~/ligesis-pcs"

SCHEMES = {
    "ligesis": {"bench_name": "ligesis_bench", "display_name": "LigeSIS"},
    "deepfold": {"bench_name": "deepfold_bench", "display_name": "DeepFold"},
    "ligero": {"bench_name": "ligero_bench", "display_name": "Ligero"},
}

DEFAULT_MUS = [24, 26, 28, 30]
DEFAULT_SCHEMES = ["ligesis", "deepfold", "ligero"]
DEFAULT_ITERATIONS = 1

# ============== 数据结构 ==============

@dataclass
class BenchResult:
    scheme: str
    mu: int
    iteration: int
    timestamp: str
    success: bool
    setup_time_ms: Optional[float] = None
    commit_time_ms: Optional[float] = None
    open_time_ms: Optional[float] = None
    verify_time_ms: Optional[float] = None
    prover_time_ms: Optional[float] = None
    total_time_ms: Optional[float] = None
    raw_output: str = ""
    error: str = ""


# ============== 服务器管理 ==============

def load_servers_config():
    with open(SERVERS_CONFIG) as f:
        return json.load(f)


def get_server_name() -> str:
    config = load_servers_config()
    return config["servers"][0]["name"]


def get_user() -> str:
    config = load_servers_config()
    return config.get("user", "")


def run_gcloud(args: list[str], timeout: int = 300) -> subprocess.CompletedProcess:
    cmd = ["gcloud"] + args
    return subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)


def gcloud_ssh(command: str, timeout: int = 3600) -> dict:
    instance = get_server_name()
    user = get_user()

    cmd = ["gcloud", "compute", "ssh"]
    if user:
        cmd.append(f"{user}@{instance}")
    else:
        cmd.append(instance)

    control_path = f"/tmp/ssh-{instance}-%r@%h:%p"
    cmd.extend([
        "--zone", ZONE,
        "--ssh-flag=-o ControlMaster=auto",
        f"--ssh-flag=-o ControlPath={control_path}",
        "--ssh-flag=-o ControlPersist=600",
        "--", "-T", command
    ])

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        return {"stdout": result.stdout, "stderr": result.stderr, "returncode": result.returncode}
    except subprocess.TimeoutExpired:
        return {"stdout": "", "stderr": "Timeout", "returncode": -1}
    except Exception as e:
        return {"stdout": "", "stderr": str(e), "returncode": -1}


def gcloud_scp(local_path: str, remote_path: str) -> bool:
    instance = get_server_name()
    user = get_user()

    cmd = ["gcloud", "compute", "scp", "--zone", ZONE]
    remote_spec = f"{user}@{instance}:{remote_path}" if user else f"{instance}:{remote_path}"
    cmd.extend([local_path, remote_spec])

    result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    return result.returncode == 0


def cmd_status(_args=None):
    result = run_gcloud(["compute", "instances", "list", "--filter=name~'^node-'"])
    print(result.stdout)
    return 0


def cmd_start(_args=None):
    server = get_server_name()
    print(f"启动服务器: {server}")

    result = run_gcloud(["compute", "instances", "start", server, "--zone", ZONE], timeout=180)
    print(result.stdout)
    if result.returncode != 0:
        print(f"启动失败: {result.stderr}", file=sys.stderr)
        return 1

    print("等待服务器就绪...")
    time.sleep(15)
    print("服务器已启动")
    return 0


def cmd_stop(_args=None):
    result = run_gcloud(["compute", "instances", "list",
                         "--filter=name~'^node-' AND status=RUNNING",
                         "--format=value(name)"])
    running = [s for s in result.stdout.strip().split('\n') if s]

    if not running:
        print("没有运行中的服务器")
        return 0

    print(f"停止 {len(running)} 个服务器: {', '.join(running)}")
    run_gcloud(["compute", "instances", "stop"] + running + ["--zone", ZONE], timeout=180)
    print("服务器已停止")
    return 0


def cmd_sync(_args=None):
    server = get_server_name()
    print(f"同步代码到 {server}...")

    tar_path = f"/tmp/ligesis_sync_{datetime.now().strftime('%H%M%S')}.tar.gz"
    exclude_args = ["--exclude=target", "--exclude=.git", "--exclude=bench_results"]

    tar_cmd = ["tar", "czf", tar_path] + exclude_args + ["-C", str(WORKSPACE), "."]
    result = subprocess.run(tar_cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"创建 tarball 失败: {result.stderr}")
        return 1

    print("  上传中...")
    if not gcloud_scp(tar_path, "~/ligesis_sync.tar.gz"):
        print("上传失败")
        return 1

    print("  解压中...")
    result = gcloud_ssh(
        f"mkdir -p {REMOTE_DIR} && cd {REMOTE_DIR} && rm -rf * && "
        f"tar xzf ~/ligesis_sync.tar.gz && rm ~/ligesis_sync.tar.gz",
        timeout=120
    )
    if result["returncode"] != 0:
        print(f"解压失败: {result['stderr']}")
        return 1

    Path(tar_path).unlink(missing_ok=True)
    print("✓ 同步完成")
    return 0


# ============== 解析器 ==============

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
    patterns = {
        "setup": [r'Setup[:\s]+([^\n]+)'],
        "commit": [r'Commit[:\s]+([^\n]+)', r'Commit \(x\d+\)[:\s]+[^\(]+\(avg[:\s]+([^\)]+)\)'],
        "open": [r'Open[:\s]+([^\n]+)', r'Open \(x\d+\)[:\s]+[^\(]+\(avg[:\s]+([^\)]+)\)'],
        "verify": [r'Verify[:\s]+([^\n]+)', r'Verify \(x\d+\)[:\s]+[^\(]+\(avg[:\s]+([^\)]+)\)'],
        "total": [r'Total[:\s]+([^\n]+)'],
    }

    for key, pats in patterns.items():
        for pat in pats:
            m = re.search(pat, output, re.IGNORECASE)
            if m:
                duration = parse_duration(m.group(1))
                if duration is not None:
                    result[key] = duration
                    break
    return result


# ============== 运行 Benchmark ==============

def check_bench_exists(scheme: str) -> bool:
    bench_name = SCHEMES[scheme]["bench_name"]
    bench_path = WORKSPACE / "ligesis-pcs" / "benches" / f"{bench_name}.rs"
    return bench_path.exists()


def run_benchmark_remote(scheme: str, mu: int, iterations: int = 1) -> BenchResult:
    config = SCHEMES.get(scheme)
    if not config:
        return BenchResult(
            scheme=scheme, mu=mu, iteration=1,
            timestamp=datetime.now().isoformat(),
            success=False, error=f"未知方案: {scheme}"
        )

    bench_name = config["bench_name"]
    cmd = (
        f"cd {REMOTE_DIR} && source ~/.cargo/env && "
        f"cargo bench --package ligesis-pcs --bench {bench_name} "
        f"--features print-trace -- --mu {mu} --iterations {iterations} 2>&1"
    )

    print(f"  执行中...", end="", flush=True)
    result = gcloud_ssh(cmd, timeout=3600)
    output = result["stdout"] + result["stderr"]

    if result["returncode"] != 0:
        print(" 失败")
        return BenchResult(
            scheme=scheme, mu=mu, iteration=iterations,
            timestamp=datetime.now().isoformat(),
            success=False, error=output[-500:], raw_output=output
        )

    parsed = parse_benchmark_output(output)
    if not parsed:
        print(" 解析失败")
        return BenchResult(
            scheme=scheme, mu=mu, iteration=iterations,
            timestamp=datetime.now().isoformat(),
            success=False, error="无法解析输出", raw_output=output
        )

    commit_ms = parsed.get("commit")
    open_ms = parsed.get("open")
    prover_ms = (commit_ms or 0) + (open_ms or 0) if commit_ms or open_ms else None

    print(" 完成")
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


# ============== 结果管理 ==============

def save_result(result: BenchResult, batch_id: Optional[str] = None):
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    if batch_id:
        result_file = RESULTS_DIR / f"batch_{batch_id}.jsonl"
    else:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        result_file = RESULTS_DIR / f"{result.scheme}_mu{result.mu}_{timestamp}.json"

    with open(result_file, 'a') as f:
        f.write(json.dumps(asdict(result)) + "\n")
    return result_file


def load_all_results() -> list[BenchResult]:
    results = []
    if not RESULTS_DIR.exists():
        return results

    for f in RESULTS_DIR.glob("*.json*"):
        with open(f) as fp:
            for line in fp:
                line = line.strip()
                if line:
                    try:
                        data = json.loads(line)
                        results.append(BenchResult(**data))
                    except json.JSONDecodeError:
                        pass
    return results


def format_time(ms: Optional[float]) -> str:
    if ms is None:
        return "-"
    if ms >= 1000:
        return f"{ms/1000:.2f}s"
    return f"{ms:.1f}ms"


# ============== 命令处理 ==============

def cmd_run(args):
    scheme = args.scheme.lower()
    mu = args.mu
    iterations = args.iterations

    print(f"\n{'='*50}")
    print(f"Running: {SCHEMES[scheme]['display_name']}, mu={mu}")
    print(f"{'='*50}\n")

    result = run_benchmark_remote(scheme, mu, iterations)

    if result.success:
        print(f"\n✓ 成功")
        print(f"  Setup:   {format_time(result.setup_time_ms)}")
        print(f"  Commit:  {format_time(result.commit_time_ms)}")
        print(f"  Open:    {format_time(result.open_time_ms)}")
        print(f"  Verify:  {format_time(result.verify_time_ms)}")
        print(f"  Prover:  {format_time(result.prover_time_ms)}")
    else:
        print(f"\n✗ 失败: {result.error[:200]}")

    result_file = save_result(result)
    print(f"\n结果已保存到: {result_file}")
    return 0 if result.success else 1


def cmd_batch(args):
    schemes = [s.strip().lower() for s in args.schemes.split(',')] if args.schemes else DEFAULT_SCHEMES
    mus = [int(m.strip()) for m in args.mus.split(',')] if args.mus else DEFAULT_MUS
    iterations = args.iterations

    available_schemes = []
    for scheme in schemes:
        if scheme not in SCHEMES:
            print(f"⚠ 未知方案: {scheme}, 跳过")
            continue
        if not check_bench_exists(scheme):
            print(f"⚠ {scheme} benchmark 不存在, 跳过")
            continue
        available_schemes.append(scheme)

    if not available_schemes:
        print("没有可用的 benchmark")
        return 1

    tests = [(scheme, mu) for scheme in available_schemes for mu in mus]
    total = len(tests)

    print(f"\n{'='*60}")
    print(f"Single-Thread PCS Benchmark (Remote)")
    print(f"服务器: {get_server_name()}")
    print(f"方案: {', '.join(available_schemes)}")
    print(f"mu: {', '.join(map(str, mus))}")
    print(f"迭代次数: {iterations}")
    print(f"总测试数: {total}")
    print(f"{'='*60}\n")

    batch_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    results = []

    for i, (scheme, mu) in enumerate(tests, 1):
        print(f"[{i}/{total}] {SCHEMES[scheme]['display_name']} mu={mu}")
        result = run_benchmark_remote(scheme, mu, iterations)
        results.append(result)
        save_result(result, batch_id)

        if result.success:
            print(f"        prover={format_time(result.prover_time_ms)}, verify={format_time(result.verify_time_ms)}\n")
        else:
            print(f"        错误: {result.error[:80]}\n")

    print(f"\n{'='*60}")
    print("测试完成! 结果汇总:")
    print(f"{'='*60}\n")
    print_summary_table(results)

    result_file = RESULTS_DIR / f"batch_{batch_id}.jsonl"
    print(f"\n结果已保存到: {result_file}")
    return 0


def cmd_report(args):
    results = load_all_results()
    if not results:
        print("没有找到结果数据")
        return 1

    successful = [r for r in results if r.success]
    if not successful:
        print("没有成功的测试结果")
        return 1

    latest = {}
    for r in sorted(successful, key=lambda x: x.timestamp):
        latest[(r.scheme, r.mu)] = r

    results_list = sorted(latest.values(), key=lambda x: (x.scheme, x.mu))
    print_summary_table(results_list)

    if args.csv:
        export_csv(results_list, args.csv)
        print(f"\nCSV 已导出到: {args.csv}")
    return 0


def print_summary_table(results: list[BenchResult]):
    by_scheme = {}
    for r in results:
        by_scheme.setdefault(r.scheme, []).append(r)

    all_mus = sorted(set(r.mu for r in results))

    header = f"{'Scheme':<12} {'mu':>4} | {'Commit':>10} {'Open':>10} {'Prover':>10} | {'Verify':>10}"
    print(header)
    print("-" * len(header))

    for scheme in sorted(by_scheme.keys()):
        scheme_results = {r.mu: r for r in by_scheme[scheme]}
        display_name = SCHEMES.get(scheme, {}).get("display_name", scheme)
        for mu in all_mus:
            r = scheme_results.get(mu)
            if r:
                print(f"{display_name:<12} {mu:>4} | "
                      f"{format_time(r.commit_time_ms):>10} "
                      f"{format_time(r.open_time_ms):>10} "
                      f"{format_time(r.prover_time_ms):>10} | "
                      f"{format_time(r.verify_time_ms):>10}")


def export_csv(results: list[BenchResult], filename: str):
    import csv
    with open(filename, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Scheme", "mu", "Commit (ms)", "Open (ms)", "Prover (ms)", "Verify (ms)", "Total (ms)", "Timestamp"])
        for r in results:
            writer.writerow([
                SCHEMES.get(r.scheme, {}).get("display_name", r.scheme),
                r.mu,
                f"{r.commit_time_ms:.2f}" if r.commit_time_ms else "",
                f"{r.open_time_ms:.2f}" if r.open_time_ms else "",
                f"{r.prover_time_ms:.2f}" if r.prover_time_ms else "",
                f"{r.verify_time_ms:.2f}" if r.verify_time_ms else "",
                f"{r.total_time_ms:.2f}" if r.total_time_ms else "",
                r.timestamp,
            ])


def cmd_list(_args=None):
    print("\n可用的 Benchmark:")
    print("-" * 40)
    for scheme, config in SCHEMES.items():
        exists = "✓" if check_bench_exists(scheme) else "✗"
        print(f"  {exists} {config['display_name']:<12} ({config['bench_name']})")
    print()
    return 0


def cmd_help(_args=None):
    print("""
可用命令:
  status              查看服务器状态
  start               启动服务器
  stop                停止服务器
  sync                同步代码到服务器
  list                列出可用 benchmark

  run -s <scheme> -m <mu> [-i <iterations>]
                      运行单个测试
                      示例: run -s ligesis -m 24

  batch [-s <schemes>] [-m <mus>] [-i <iterations>]
                      批量运行测试
                      示例: batch -s ligesis,deepfold -m 24,26,28

  report [--csv <file>]
                      查看/导出结果

  help                显示此帮助
  exit, quit, q       退出
""")
    return 0


# ============== 交互式模式 ==============

def create_parser():
    parser = argparse.ArgumentParser(prog="", add_help=False)
    subparsers = parser.add_subparsers(dest="command")

    subparsers.add_parser("status")
    subparsers.add_parser("start")
    subparsers.add_parser("stop")
    subparsers.add_parser("sync")
    subparsers.add_parser("list")
    subparsers.add_parser("help")

    p = subparsers.add_parser("run")
    p.add_argument("--scheme", "-s", type=str, required=True, choices=list(SCHEMES.keys()))
    p.add_argument("--mu", "-m", type=int, default=24)
    p.add_argument("--iterations", "-i", type=int, default=1)

    p = subparsers.add_parser("batch")
    p.add_argument("--schemes", "-s", type=str, default=None)
    p.add_argument("--mus", "-m", type=str, default=None)
    p.add_argument("--iterations", "-i", type=int, default=DEFAULT_ITERATIONS)

    p = subparsers.add_parser("report")
    p.add_argument("--csv", type=str, default=None)

    return parser


def interactive_mode():
    print("=" * 60)
    print("Single-Thread PCS Benchmark - 交互式模式")
    print("输入 'help' 查看可用命令, 'exit' 退出")
    print("=" * 60)

    parser = create_parser()
    cmd_map = {
        "status": cmd_status,
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
            print("\n退出")
            break

        if not line:
            continue

        if line.lower() in ("exit", "quit", "q"):
            print("退出")
            break

        try:
            argv = shlex.split(line)
            args = parser.parse_args(argv)

            if args.command in cmd_map:
                cmd_map[args.command](args)
            else:
                print(f"未知命令: {line}")
                cmd_help()

        except SystemExit:
            # argparse 在解析失败时会调用 sys.exit()
            pass
        except Exception as e:
            print(f"错误: {e}")


# ============== 主程序 ==============

def main():
    # 如果没有参数，进入交互式模式
    if len(sys.argv) == 1:
        interactive_mode()
        return 0

    # 否则执行单次命令
    parser = argparse.ArgumentParser(
        description="Single-Thread PCS Benchmark (Remote Server)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  %(prog)s                                         # 交互式模式
  %(prog)s status                                  # 查看服务器状态
  %(prog)s start                                   # 启动服务器
  %(prog)s sync                                    # 同步代码
  %(prog)s run --scheme ligesis --mu 24            # 运行单个测试
  %(prog)s batch --schemes ligesis,deepfold --mus 24,26
  %(prog)s report --csv results.csv                # 导出结果
        """
    )

    subparsers = parser.add_subparsers(dest="command", help="命令")

    subparsers.add_parser("status", help="查看服务器状态").set_defaults(func=cmd_status)
    subparsers.add_parser("start", help="启动服务器").set_defaults(func=cmd_start)
    subparsers.add_parser("stop", help="停止服务器").set_defaults(func=cmd_stop)
    subparsers.add_parser("sync", help="同步代码到服务器").set_defaults(func=cmd_sync)
    subparsers.add_parser("list", help="列出可用的 benchmark").set_defaults(func=cmd_list)

    p = subparsers.add_parser("run", help="运行单个测试")
    p.add_argument("--scheme", "-s", type=str, required=True, choices=list(SCHEMES.keys()))
    p.add_argument("--mu", "-m", type=int, default=24)
    p.add_argument("--iterations", "-i", type=int, default=1)
    p.set_defaults(func=cmd_run)

    p = subparsers.add_parser("batch", help="批量运行测试")
    p.add_argument("--schemes", "-s", type=str, default=None,
                   help=f"方案列表，逗号分隔 (默认: {','.join(DEFAULT_SCHEMES)})")
    p.add_argument("--mus", "-m", type=str, default=None,
                   help=f"mu 列表，逗号分隔 (默认: {','.join(map(str, DEFAULT_MUS))})")
    p.add_argument("--iterations", "-i", type=int, default=DEFAULT_ITERATIONS)
    p.set_defaults(func=cmd_batch)

    p = subparsers.add_parser("report", help="查看/导出结果")
    p.add_argument("--csv", type=str, default=None)
    p.set_defaults(func=cmd_report)

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        return 1

    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
