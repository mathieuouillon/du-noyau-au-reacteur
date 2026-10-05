"""precalcul_micro.py -- Calcule une fois pour toutes l'energie microscopique
E_micro(delta) = dE_couches + dE_appariement de chaque noyau, sur la grille de
deformations. Elle ne depend d'aucun parametre ajuste : on la met en cache
dans micro.npz, que etude_frdm.py relit.

Lancer :  python precalcul_micro.py            (tout, ~15 min sur un coeur)
          python precalcul_micro.py 0 500      (une tranche, reprise possible)
"""
import sys, os, time, glob
import numpy as np
from donnees_liaison import table_liaison
from frdm import TablesMicro, GRILLE_DELTA


def liste_noyaux():
    t = table_liaison()                     # mesures ET estimations
    t = t[(t["Z"] >= 8) & (t["N"] >= 8)]
    return list(zip(t["Z"].tolist(), t["N"].tolist()))


def tranche(debut, fin):
    zn = liste_noyaux()[debut:fin]
    T = TablesMicro()
    avec = [[T.energie_micro(z, n, d, True) for d in GRILLE_DELTA] for z, n in zn]
    sans = [[T.energie_micro(z, n, d, False) for d in GRILLE_DELTA] for z, n in zn]
    np.savez(f"micro_{debut:05d}.npz", Z=[a for a, _ in zn], N=[b for _, b in zn],
             micro=avec, couches_seules=sans)


def assembler():
    parts = sorted(glob.glob("micro_0*.npz"))
    d = [np.load(p) for p in parts]
    np.savez("micro.npz", delta=GRILLE_DELTA,
             Z=np.concatenate([x["Z"] for x in d]), N=np.concatenate([x["N"] for x in d]),
             micro=np.concatenate([x["micro"] for x in d]),
             couches_seules=np.concatenate([x["couches_seules"] for x in d]))
    for p in parts:
        os.remove(p)


if __name__ == "__main__":
    if len(sys.argv) == 3:
        tranche(int(sys.argv[1]), int(sys.argv[2]))
    else:
        n = len(liste_noyaux())
        for k in range(0, n, 300):
            t0 = time.time()
            tranche(k, k + 300)
            print(f"  {min(k + 300, n)}/{n}  ({time.time() - t0:.0f} s)")
        assembler()
        print("micro.npz ecrit")
