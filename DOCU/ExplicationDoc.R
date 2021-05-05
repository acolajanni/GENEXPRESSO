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

############################################################################################################################################################################

# Package Renv explication

# install.packages("renv")

# Objectif : créer une libraire dite locale, pour chaque projet. Ex : projet 1 / libraire 1 ; projet 2 / librairie 2. On "décentralise" la bibliothèque en somme.

# Ca permet de changer la librairie d'un projet sans impacter les autres projets (= isolation), ensuite la sauvegarde fixe de vos packages sont transmissibles aux autres personnes
# voulant bosser sur le projet (= portabilité), et enfin on peut réutiliser la sauvegarde des packages (sauvegardée dans un lockfile) et restaurer la version sauvegardée au cas où, un peu comme quand j'ai fait failli
# tout faire capoter sur Git l'autre fois et que vous avez pu restaurer la dernière version (= reproductibilité).

# Etape 1 : initialiser l'environnement local de votre projet
# renv::init()
# A partir de maintenant, tout packages installés et chargés se feront dans cette bibliothèque locale, donc biblio 1 du projet 1 par exemple.
# Pour vérifier, faites   .libPaths()  et regarder si dans les 2 chemins en sortie il y a bien /renv/ dans le chemin de biblio "perso [1] et /renv-system-library dans le [2]

# Etape 2 : sauvegarder
# Vous avez installé de nouveau package dans votre librairie 1 du projet 1. Maintenant il faut sauvegarder cette modification :
# renv::snapshot()
# Il va vous dire quels packages vont être ajoutés dans la sauvegarde fixe et leur version, il faudra ensuite confirmer dans la console.
# Il va stocker le nouvel état de la biblio dans un fichier nommé renv.lock, contenant les packages, leurs sources etc
# Si j'ai bien capté c'est ce qu'attends les profs, car c'est avec ça qu'on va pouvoir transmettre l'état des packages fonctionnels pour le projet (la portabilité).

# Etape 3 : restaurer
# Mauvaises installations ou autres modifications sur les packages non désirées ? pas de panique, on peut remettre l'état de la biblio avec la dernière sauvegarde faite avec ::snap'
# renv::restore()
# Ca va récupérer le fichier lock pour l'appliquer à notre biblio, retour arrière. 
