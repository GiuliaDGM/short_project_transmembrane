# SUJET : ASSIGNATION ET DÉTECTION DES PARTIES TRANSMEMBRANAIRES D'UNE PROTÉINE

**Objectif** : Réaliser un programme reprenant la méthode décrite dans l'article mentionné afin de concevoir un outil permettant de déterminer les zones transmembranaires d'une protéine.

### Étapes :

1. Choisir des protéines membranaires et globulaires de référence (simples à analyser, voir la banque OPM).
2. Calculer la surface accessible au solvant de chaque acide aminé à l'aide de **DSSP**.
3. Lire le fichier PDB et extraire les coordonnées des atomes **C-alpha**.
4. Calculer le centre de masse de la protéine.
5. Déterminer les droites passant par le centre de masse et quadriller avec suffisamment de précision toutes les directions.
6. Déplacer une tranche de 1 Angström normale à une droite et calculer l'hydrophobicité relative des résidus exposés dans la tranche.
7. Calculer la position de la membrane en faisant la moyenne de l'hydrophobicité relative et en comparant ces valeurs selon les différentes droites.

**Référence** : Tusnády GE, Dosztányi Z, Simon I. *Transmembrane proteins in the Protein Data Bank: identification and classification*. Bioinformatics. 2004 Nov 22;20(17):2964-72. Epub 2004 Jun 4.

**Remarques** :

- Il n'est pas demandé de vérifier si la protéine est de taille suffisante ou si le fichier contient effectivement une protéine (on prendra donc soin de vérifier son programme sur des chaînes uniques de protéines membranaires et de protéines globulaires).
- Il n'est pas demandé de reconstruire le modèle biologique de la molécule (à partir du champ BIOMOLECULE).
- Le programme **DSSP** sera utilisé pour calculer la surface accessible au solvant.
