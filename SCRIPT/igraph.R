d1 <- rbind(c(" A ", "B "),c(" A ", "C "), c(" B ", "C "), c(" A ", "D "), c(" D ", "E "))
d2 <- rbind(c(" A ", "F "), c(" A ", "D "), c(" A ", "C "), c(" C ", "I "))

g1 <- graph.data.frame(d1, directed = FALSE)
g2 <- graph.data.frame(d2, directed = FALSE)

#g1 = spearman
#g2 = kendall


par(mfrow = c(1,2))
plot(g1, main =  "g1" )
plot(g2, main = "g2")

#3 points communs (ACD) et deux liens communs (AD et AC)

#ajout d’attributs : degré des sommets
V(g1)$degree <-  degree(g1)
V(g2)$degree <- degree(g2)

#ajout d’attributs aux liens (intermédiarité)
E(g1)$between <- edge.betweenness(g1)
E(g2)$between <- edge.betweenness(g2)

#ajout d’attributs au réseau : densité
g1$densite <- graph.density(g1)
g2$densite <- graph.density(g2)

#fonction graph.complementer
#crée un graphe avec les liens absents du graphe de départ
#les liens présents au départ sont eux supprimés
compg1 <- graph.complementer(g1, loops = FALSE)
plot(compg1)
#les attributs des sommets et du graphe sont conservés
#attributs des liens sont perdus
V(compg1)$degree


#fonction graph.difference premier réseau – deuxième réseau (opérateur %m%)
#seuls les liens présents dans le premier réseau et non dans le deuxième sont conservés
#tous les attributs du premier réseau sont gardés
diffg1 <- graph.difference(g1, g2)
diffg1

diffg2 <- graph.difference(g2, g1)

#coordonnées (quasi) identiques pour faciliter la comparaison
lay <- layout.fruchterman.reingold(g1)

#lien A-D présent dans g1 et g2 donc g1-g2 et g2-g1 entrainent sa disparition
par(mfrow = c(2,2))
plot(g1, main = "g1", layout = lay, vertex.color = "red", vertex.size = 20, vertex.label.dist =3)
plot(g2, main = "g2", layout = lay, vertex.color = "red", vertex.size = 20, vertex.label.dist =3)
plot(diffg1, main = "graph.difference(g1, g2)", layout = lay, vertex.color = "red", vertex.size = 20, vertex.label.dist =3)
plot(diffg2, main = "graph.difference(g2, g1)", layout = lay, vertex.color = "red", vertex.size = 20, vertex.label.dist =3)
dev.off()

#calculer degré actualisé
V(diffg1)$degree2 <- degree(diffg1)
V(diffg1)$degree
V(diffg1)$degree2

#graph.disjoint.union pour joindre des graphes dont les sommets différents (opérateur %du%)

#liens communs à deux graphes : conserve les attributs
#tous les sommets sont conservés (opérateur %s%)
interg1 <- graph.intersection(g1,g2, keep.all.vertices = TRUE)

#seuls sommets adjacents aux liens communs sont gardés (ici ACD)
interg2 <- graph.intersection(g1,g2, keep.all.vertices = FALSE)

par(mfrow = c(1,2))
plot(interg1, sub = "graph.intersection(g1,g2,\nkeep.all.vertices = TRUE)")
plot(interg2, sub = "graph.intersection(g1,g2,\nkeep.all.vertices = FALSE)")
dev.off()

#union de graphes : tout lien présent dans un moins un graphe
#est dans le graphe produit – tous les attributs sont conservés
#ne crée pas de liens multiples

uniong <- graph.union(g1,g2) #opérateur %u%
plot(uniong)



## A simple example with a couple of actors
## The typical case is that these tables are read in from files....
actors <- data.frame(name=c("Alice", "Bob", "Cecil", "David",
                            "Esmeralda"))
relations <- data.frame(from=c("Bob", "Cecil", "Cecil", "David",
                               "David", "Esmeralda"),
                        to=c("Alice", "Bob", "Alice", "Alice", "Bob", "Alice"))
g <- graph_from_data_frame(relations, directed=F, vertices=actors)
print(g, e=TRUE, v=TRUE)
plot(g)
## The opposite operation
as_data_frame(g, what="vertices")
as_data_frame(g, what="edges")






