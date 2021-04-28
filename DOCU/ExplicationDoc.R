#Create doc

####### 1 ERE VERSION (possibilité que j'affine avec mes recherches, mais ce qui a ici devrait faire l'affaire)

# New project => Create R package => saisir le nom du package
# Ouvrir le script avec les fonctions
# Build => Configure build tools => generate document with roxygen
# Aller dans le fichier DESCRIPTIONS, mettre le titre du package, les auteurs,...licence ( = GPL-3 ou MIT). Save + fermer
# Ouvrir votre script, puis aller DANS une fonction => Code => Insert Roxygen skeleton
# Mettre le titre, les paramètres, ce que ça return, avec pourquoi pas un exemple d'utilisation
# Puis faire Install et restart (Build) => ensuite faire More => Document
# Fichier man : liste des docs des fonctions du package créés (faire Preview pour voir une doc)
# Pour voir la tête du package avec la doc de novuelles fonctions (fichier man), faire Install et restart (il permettra d'être utilisé de manière générale dans R)
# Puis aller dans Help => Home => Packages => Trouver le nom du package et cliquer dessus.
# FIN

# Essayer avec library()
# NOTE : on peut utiliser Check dans Build pour vérifier si tout va bien dans la construction de notre package
