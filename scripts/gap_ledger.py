#!/usr/bin/env python3
"""gap_ledger.py <missing_described.txt> -- tally the giac gap ledger's buckets (crates/cas-giac/README.md ## Tests).

Buckets of missing_described.txt (gap_describe.py missing), from R1.8-003, re-checked in 1.8 T2:
N new module, P plotting, S syntax, I I/O / TI screens / GUI, M ML/music/geodesy,
D already a README departure, A alias/inert spelling/option of a built name, C computational, existing module.
C, A, D and M are keyed by an entry's first name; N, P, S and I by line ranges of the 2026-11-26 list
(giac 1.9.0 aide_cas, en.txt at 1.8 T2). The script fails when a line has no bucket or two, or a named entry
is missing, so a changed list (new giac, new names) shows up instead of shifting the counts.
"""
import sys

NAMED = {
    "C": "REDIM back taux_accroissement shift_phase linabs GF stft istft hht envelope border non_recursive_normal "
         "semi_augment root range roots simult newMat Fill avgRC multinomial shift upper lower format additionally "
         "atan2 a2q basis changebase has lin moyal order q2a choice hasard suppress trn truncate divis fmod modf "
         "dayofweek period boxcar tri sinc addtable box_constraints",
    "A": "besselY besselJ fadeev rdiv maple_ifactors Int subs getType list polynom parfrac base array Inverse at "
         "Det ' Hilbert /% %/ div Nullspace expln ! Rref | mod Fourier",
    "D": "QR LU confrac Factor Gcd Resultant xcas_mode python_compat Quo Rem gbasis_reinject rur_gbasis",
    "M": "tonnetz neural_network train MSE geodetic2ecef ecef2geodetic",
}
RANGES = {
    "N": "13-16 31 42-43 50 53-68 70-73 107-127 131-341 346-350 352 355 415-417 448 454 464-470 479-480 485 490 "
         "493 501-502 510-511 515-519 531 534 537-541 550-555 557-559 563-568 573 577 599-600 607 615-616 619-621 "
         "650 658-659 661 667-672 704 707 710-714 763-764 774-776 779-975 977-997",
    "P": "3 7-8 17 20-29 39-40 44 48 51-52 69 74 82-105 342-345 353 391-395 403-405 408-410 433 436 455-459 "
         "472-474 481-484 486-489 491-492 535 579-596 618 622-627 629-647 666 695-703 731-746 752 762 767-768",
    "S": "4 18-19 32 37 41 75-81 128-130 360 369-370 378-388 390 396 406 412 430 447 450-451 453 494 499 503-505 "
         "507-508 513 521-522 524 526 528-530 532-533 536 543-544 560-562 598 608 610 612-614 617 651-652 665 "
         "673 676-691 694 709 717-719 721-723 727-730 747-748",
    "I": "30 33-36 45 351 356-358 361-366 371-376 389 397-402 420-421 425-427 429 434 437-440 442-446 452 506 509 "
         "514 520 523 546-548 576 603-606 648-649 653 662-663 674-675 705-706 708 720 753-761 777-778 1004",
}


def lines_of(spec):
    for part in spec.split():
        a, _, b = part.partition("-")
        yield from range(int(a), int(b or a) + 1)


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    rows = [l.rstrip("\n") for l in open(sys.argv[1], encoding="utf-8")]
    first = [r.split(" | ")[0].split()[0] if r.split(" | ")[0].split() else r.split(" | ")[0] for r in rows]
    where = {}
    for b, spec in RANGES.items():
        for n in lines_of(spec):
            where.setdefault(n, []).append(b)
    bad = 0
    for b, names in NAMED.items():
        for name in names.split():
            hits = [i + 1 for i, f in enumerate(first) if f == name]
            if not hits:
                print(f"named {b} entry {name!r} not in the list"); bad += 1
            for n in hits:
                where.setdefault(n, []).append(b)
    for n in range(1, len(rows) + 1):
        bs = where.get(n, [])
        if len(bs) != 1:
            print(f"line {n} ({first[n - 1]}): buckets {bs or 'none'}"); bad += 1
    extra = [n for n in where if n > len(rows)]
    if extra:
        print(f"ranges past the end: {extra[:5]}"); bad += 1
    counts = {}
    for n in range(1, len(rows) + 1):
        for b in where.get(n, []):
            counts[b] = counts.get(b, 0) + 1
    print(" ".join(f"{b}={counts.get(b, 0)}" for b in "NPSIMDAC"), f"total={len(rows)}")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
