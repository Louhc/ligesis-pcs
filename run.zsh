# python3 benchmark.py sync

# gtimeout 3000s python3 benchmark.py batch -s dbasefold-pcs -m 23,24,25,26 -i 10 --build
# gtimeout 3000s python3 benchmark.py batch -s dhyperfond -m 19,20,21,22,23,24,25 -i 10 --build
# gtimeout 3000s python3 benchmark.py batch -s dhyperpianist -m 19,20,21,22,23 -i 10 --build
# gtimeout 3000s python3 benchmark.py batch -s dsumcheck3 -m 26,27,28,29,30 -i 10 --build
# gtimeout 3000s python3 benchmark.py batch -s dsumcheck4 -m 26,27,28,29,30 -i 10 --build

# python3 benchmark.py stop
# python3 benchmark.py set-n 4
# python3 benchmark.py set-vm 1-4 64g
# python3 benchmark.py start
# gtimeout 3000s python3 benchmark.py batch -s dbasefold-pcs -m 23,24,25,26 -i 10 --build
# gtimeout 3000s python3 benchmark.py batch -s ddory-pcs -m 11,12,13,14,15,16,17,18 -i 10 --build
# gtimeout 3000s python3 benchmark.py batch -s dhyperfond -m 18,19,20,21,22,23,24,25 -i 10 --build
# gtimeout 3000s python3 benchmark.py batch -s dhyperpianist -m 18,19,20,21,22,23 -i 10 --build
# gtimeout 3000s python3 benchmark.py batch -s dmkzg-pcs -m 19,20,21,22,23,24 -i 10 --build
# gtimeout 3000s python3 benchmark.py batch -s dpip-fri-pcs -m 24,25,26,27,28 -i 10 --build
# gtimeout 3000s python3 benchmark.py batch -s dsumcheck3 -m 26,27,28,29,30 -i 10 --build
# gtimeout 3000s python3 benchmark.py batch -s dsumcheck3 -m 27,29 -i 10 --build
# gtimeout 3000s python3 benchmark.py batch -s dsumcheck4 -m 26,27,28,29 -i 10 --build

python3 benchmark.py stop
python3 benchmark.py set-n 2
# python3 benchmark.py set-vm 1-4 64g
python3 benchmark.py start
gtimeout 3000s python3 benchmark.py batch -s dbasefold-pcs -m 22,23,24,25,26 -i 10 --build
gtimeout 3000s python3 benchmark.py batch -s ddory-pcs -m 10,11,12,13,14,15,16,17,18 -i 10 --build
gtimeout 3000s python3 benchmark.py batch -s dhyperfond -m 17,18,19,20,21,22,23,24,25 -i 10 --build
gtimeout 3000s python3 benchmark.py batch -s dhyperpianist -m 17,18,19,20,21,22,23 -i 10 --build
gtimeout 3000s python3 benchmark.py batch -s dmkzg-pcs -m 18,19,20,21,22,23,24 -i 10 --build
gtimeout 3000s python3 benchmark.py batch -s dpip-fri-pcs -m 23,24,25,26,27,28 -i 10 --build
gtimeout 3000s python3 benchmark.py batch -s dsumcheck3 -m 26,27,28,29,30 -i 10 --build
gtimeout 3000s python3 benchmark.py batch -s dsumcheck4 -m 26,27,28,29,30 -i 10 --build

python3 benchmark.py stop