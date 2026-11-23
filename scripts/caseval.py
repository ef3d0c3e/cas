#!/usr/bin/env python3
"""caseval.py <inputs.cas> [timeout_s] -- giac's answer to each input line through libgiac's caseval (C0018).

Run from anywhere. Each line runs in its own child process, so a hang or crash costs one line.
Lines starting with '#' or blank are skipped. Output: `input  =>  answer`, one per line.
"""
import subprocess
import sys


def one(text):
    from ctypes import CDLL, c_char_p
    g = CDLL("libgiac.so.0")
    g.caseval.restype = c_char_p
    g.caseval.argtypes = [c_char_p]
    out = g.caseval(text.encode())
    print(out.decode(errors="replace") if out is not None else "<null>")


def main():
    if sys.argv[1] == "--one":
        one(sys.argv[2])
        return
    path = sys.argv[1]
    limit = float(sys.argv[2]) if len(sys.argv) > 2 else 30.0
    with open(path, encoding="utf-8") as f:
        lines = [l.rstrip("\n") for l in f]
    for line in lines:
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        try:
            r = subprocess.run([sys.executable, "-I", __file__, "--one", line],
                               capture_output=True, text=True, timeout=limit)
            ans = (r.stdout.strip() or ("<stderr> " + r.stderr.strip()[-200:])).replace("\n", " | ")
        except subprocess.TimeoutExpired:
            ans = f"<timeout {limit:.0f}s>"
        print(f"{line}  =>  {ans}", flush=True)


if __name__ == "__main__":
    main()
