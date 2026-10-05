"""
donnees_liaison.py -- Energies de liaison EXPERIMENTALES par nucleon.

fission.py calcule B/A avec la formule semi-empirique de masse : c'est un
MODELE. Ce module fournit les valeurs MESUREES, derivees des masses
atomiques evaluees, pour pouvoir confronter les deux.

Source des masses
-----------------
AME2020 (Atomic Mass Evaluation 2020) :
    M. Wang, W.J. Huang, F.G. Kondev, G. Audi, S. Naimi,
    "The AME 2020 atomic mass evaluation (II). Tables, graphs and references",
    Chinese Physics C 45 (2021) 030003, DOI 10.1088/1674-1137/abddaf

Acces : via le paquet Python `periodictable` (v2.1.0), qui embarque le
fichier massround.mas20.txt telecharge le 2023-07-06 depuis l'AMDC de
l'AIEA (www-nds.iaea.org/amdc/ame2020/). "massround" = version a
incertitudes arrondies, publiee avec l'evaluation.

MESURE OU ESTIMATION. Dans AME, une masse suivie de '#' n'est PAS
purement experimentale : elle est estimee par extrapolation des
tendances des noyaux voisins (systematique). On conserve ce marqueur
dans la colonne `estimee`. Sur les 2825 nuclides retenus ici, 375 sont
estimes et 2450 mesures.

Abondances naturelles : CIAAW / IUPAC, via le meme paquet. Voir
ABONDANCES_COMPLEMENT pour le cas de l'uranium.

Methode
-------
    B(A,Z) = [ Z*m(1H) + N*m(n) - M_atomique(A,Z) ] * 931.494 MeV/u

On utilise la masse ATOMIQUE de l'hydrogene 1 (noyau + 1 electron) et la
masse ATOMIQUE du nuclide : les Z electrons se compensent de part et
d'autre. On neglige ainsi l'energie de liaison des electrons, qui atteint
~0.8 MeV pour l'uranium, soit 4e-4 de B. Negligeable ici, mais a savoir.
"""

import numpy as np
import periodictable as pt

# Constantes CODATA
U_EN_MEV = 931.49410242         # MeV par unite de masse atomique
MASSE_NEUTRON = 1.00866491595   # u, masse du neutron libre
MASSE_H1 = 1.00782503207        # u, masse ATOMIQUE de l'hydrogene 1


# Le fichier de `periodictable` contient bien les abondances de l'uranium
# (CIAAW, minerai de Namibie de reference, Meija 2016), mais l'analyseur du
# paquet les perd : `U[235].abundance` renvoie 0. On les reinjecte ici avec
# les valeurs exactes du fichier. C'est le seul element concerne : les
# autres elements a abondance nulle (Tc, Pm, Po, At, Rn, Fr, Ra, Ac) n'ont
# effectivement aucun isotope stable.
ABONDANCES_COMPLEMENT = {
    ("U", 234): 0.0054,
    ("U", 235): 0.7204,
    ("U", 238): 99.2742,
}


def _marqueurs_estimation():
    """{(symbole, A): True si la masse AME est estimee ('#')}."""
    import periodictable.mass as pm
    marqueurs = {}
    for ligne in pm.isotope_mass.split("\n"):
        if not ligne.strip():
            continue
        nuclide, masse = ligne.split(",")[:2]
        _, symbole, A = nuclide.split("-")
        marqueurs[(symbole, int(A))] = "#" in masse
    return marqueurs


def table_liaison(z_max=100, exiger_abondance=False,
                  mesurees_seulement=False):
    """Construit la table des energies de liaison mesurees.

    Retourne un tableau structure avec, par nuclide :
        Z, N, A, symbole, masse (u), B (MeV), B_sur_A (MeV), abondance (%)

    exiger_abondance   : ne garder que les isotopes presents dans la nature.
    mesurees_seulement : ecarter les masses estimees ('#' dans AME2020).
    """
    marqueurs = _marqueurs_estimation()
    lignes = []
    for element in pt.elements:
        Z = element.number
        if Z < 1 or Z > z_max:
            continue
        for isotope in element:
            A = isotope.isotope
            masse = getattr(isotope, "mass", None)
            if masse is None or A is None or A < 1:
                continue
            abondance = getattr(isotope, "abundance", 0.0) or 0.0
            abondance = ABONDANCES_COMPLEMENT.get((element.symbol, A), abondance)
            if exiger_abondance and abondance <= 0.0:
                continue
            N = A - Z
            if N < 0:
                continue
            estimee = marqueurs.get((element.symbol, A), False)
            if mesurees_seulement and estimee:
                continue
            B = (Z * MASSE_H1 + N * MASSE_NEUTRON - masse) * U_EN_MEV
            lignes.append((Z, N, A, element.symbol, masse, B, B / A,
                           abondance, estimee))

    dtype = [("Z", int), ("N", int), ("A", int), ("symbole", "U3"),
             ("masse_u", float), ("B_MeV", float), ("B_sur_A_MeV", float),
             ("abondance_pct", float), ("estimee", bool)]
    table = np.array(lignes, dtype=dtype)
    return np.sort(table, order=["Z", "A"])


def chercher(table, symbole, A):
    """Retourne la ligne d'un nuclide donne, ou None."""
    masque = (table["symbole"] == symbole) & (table["A"] == A)
    return table[masque][0] if masque.any() else None


def B_sur_A(table, symbole, A):
    """B/A mesure en MeV, ou None si le nuclide est absent de la table."""
    ligne = chercher(table, symbole, A)
    return None if ligne is None else float(ligne["B_sur_A_MeV"])


def energie_separation_neutron(table, symbole, A):
    """S_n mesuree = B(A,Z) - B(A-1,Z), en MeV."""
    haut, bas = chercher(table, symbole, A), chercher(table, symbole, A - 1)
    if haut is None or bas is None:
        return None
    return float(haut["B_MeV"] - bas["B_MeV"])


def Q_fission(table, lourd, frag1, frag2):
    """Q d'une fission, a partir des energies de liaison MESUREES.

    Chaque argument est un couple (symbole, A). Les neutrons libres ont une
    energie de liaison nulle et n'entrent pas dans le bilan.
    """
    a = chercher(table, *lourd)
    b = chercher(table, *frag1)
    c = chercher(table, *frag2)
    if a is None or b is None or c is None:
        return None
    return float(b["B_MeV"] + c["B_MeV"] - a["B_MeV"])


def exporter_csv(table, chemin):
    """Ecrit la table complete au format CSV."""
    entete = ("Z,N,A,symbole,masse_atomique_u,"
              "energie_liaison_MeV,B_sur_A_MeV,abondance_pct,origine")
    with open(chemin, "w", encoding="utf-8") as fichier:
        fichier.write(entete + "\n")
        for r in table:
            fichier.write(f"{r['Z']},{r['N']},{r['A']},{r['symbole']},"
                          f"{r['masse_u']:.9f},{r['B_MeV']:.4f},"
                          f"{r['B_sur_A_MeV']:.5f},{r['abondance_pct']:.6f},"
                          f"{'estimee' if r['estimee'] else 'mesuree'}\n")
    return chemin


if __name__ == "__main__":
    t = table_liaison()
    print(f"{len(t)} nuclides")
    print(exporter_csv(t, "energies_liaison.csv"))
