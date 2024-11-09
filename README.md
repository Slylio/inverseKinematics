# Inverse Kinematics - CCD Method

## Description
Ce projet implémente une méthode d'Incinération Cinématique (IK) utilisant l'algorithme CCD (Cyclic Coordinate Descent) pour contrôler une chaîne d'articulations 2D. L'objectif est d'ajuster les angles des articulations afin que l'effecteur final atteigne une cible donnée.

Le projet utilise la méthode de la jacobienne pour ajuster les angles des articulations à chaque itération et rapproche progressivement l'effecteur final de la cible.

On pourra se servir de ce mini projet pour implémenter la cinématique inverse sur un bras robot articulé en 3D ?
https://medium.com/unity3danimation/overview-of-jacobian-ik-a33939639ab2 (approche informatique qui m'a aidée)

## Structure du projet
- **main.py** : Le script principal avant l'ajout de la jacobienne. Ce fichier m'a servi de base pour mieux comprendre le principe.
- **main_modified.py** : Ma version modifiée de `main.py`, qui inclut le calcul de la jacobienne et son utilisation pour la méthode CCD. C'est la version fonctionnelle du projet.
- **res/** : Dossier contenant les résultats de l'algorithme CCD sous forme de fichiers d'image. Chaque itération est enregistrée dans ce dossier pour suivre l'évolution de la solution.

