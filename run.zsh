# python3 benchmark.py stop
# python3 benchmark.py set-n 2
# python3 benchmark.py start

# gtimeout 3000s python3 benchmark.py batch -s dsumcheck3 -m 26,27,28,29,30 -i 10 --build
# gtimeout 3000s python3 benchmark.py batch -s dsumcheck4 -m 26,27,28,29 -i 10 --build

python3 benchmark.py stop
python3 benchmark.py set-n 8
python3 benchmark.py set-vm 1-8 32g
python3 benchmark.py start

gtimeout 3000s python3 benchmark.py batch -s dhyperfond -m 23 -i 10 --build
gtimeout 3000s python3 benchmark.py batch -s dhyperpianist -m 21 -i 10 --build
gtimeout 3000s python3 benchmark.py batch -s dsumcheck4 -m 28 -i 10 --build

python3 benchmark.py stop
python3 benchmark.py set-n 4
python3 benchmark.py set-vm 1-4 64g
python3 benchmark.py start

gtimeout 3000s python3 benchmark.py batch -s dhyperpianist -m 20 -i 10 --build

python3 benchmark.py stop