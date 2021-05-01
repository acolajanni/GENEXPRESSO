d1 <- rbind(c(" A ", "B "),c(" A ", "C "), c(" B ", "C "), c(" A ", "D "), c(" D ", "E "))
d2 <- rbind(c(" A ", "F "), c(" A ", "D "), c(" A ", "C "), c(" C ", "I "))

g1 <- graph.data.frame(d1, directed = FALSE)
g2 <- graph.data.frame(d2, directed = FALSE)

g1 = spearman
g2 = TOM
g3 = kendall
lay = layout_with_lgl(g2)


dev.off()
par(mfrow = c(1,3))
plot(g1, main =  "spearman" , lay = layout_with_lgl(g1),vertex.label = NA)
plot(g2, main = "TOM", layout = lay, vertex.label = NA)
plot(g3, main = "kendall", layout = lay)

#plot(g3, main = "g3", layout = lay)

#3 points communs (ACD) et deux liens communs (AD et AC)

#ajout d’attributs : degré des sommets
#V(g1)$degree <-  degree(g1)
#V(g2)$degree <- degree(g2)
#V(g3)$degree <- degree(g3)
#ajout d’attributs aux liens (intermédiarité)
#E(g1)$between <- edge.betweenness(g1)
#E(g2)$between <- edge.betweenness(g2)
#E(g3)$between <- edge.betweenness(g3)


unionG1_G2 = graph.union(g1,g2) #opérateur %u%
plot(unionG1_G2, layout = lay, main = "union")

#ajout d’attributs au réseau : densité
#g1$densite <- graph.density(g1)
#g2$densite <- graph.density(g2)
#g3$densite <- graph.density(g3)


# Présent dans g1 mais pas dans g2
diffg1 <- graph.difference(g1, g2)
#V(diffg1)$degree <-  degree(diffg1)
#E(diffg1)$between<-  edge.betweenness(diffg1)
#diffg1$density <-  graph.density(diffg1)
#présent dans g2 mais pas dans g1
diffg2 <- graph.difference(g2, g1)
#V(diffg2)$degree <-  degree(diffg2)
#E(diffg2)$between<-  edge.betweenness(diffg2)
#diffg2$density <-  graph.density(diffg2)
# Présent dans les deux
interg1 <- graph.intersection(g1,g2, keep.all.vertices = T)
#♦V(interg1)$degree <-  degree(interg1)
#E(interg1)$between<-  edge.betweenness(interg1)
#interg1$density <-  graph.density(interg1)



G1_edges = as_data_frame(diffg1,what = "edges")
G1_edges$color = "blue"

G2_edges = as_data_frame(diffg2,what = "edges")
G2_edges$color = "darkgreen"

Common_edges = as_data_frame(interg1, what="edges")
Common_edges = Common_edges[-c(3,4)]
Common_edges$color = "red"

test = rbind(G1_edges,G2_edges, Common_edges)


g = graph.data.frame(test, directed = F)
plot(g, layout = lay,
     edge.width = 1,
     vertex.size = 2,
     vertex.label = NA,
     vertex.color = "white",
     label.color = "black",
     label.font = 2)

legend("bottomleft",
       #x=-1.5, y=-1.1,
       c("Spearman connexions","Kendall connexions", "common connexions"), 
       pch=18, 
       col=c("blue","darkgreen","red"), 
       pt.cex=0, #taille des bulles légendes 
       cex=.8, #taille de la police légende
       lty=c(1,1,1),
       lwd = 3,
       bty="n", #absence de cadre autour de la légende 
       ncol=1)




#################################################################


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
diffg3 <- graph.difference(g1, g3)

#coordonnées (quasi) identiques pour faciliter la comparaison
lay <- layout.fruchterman.reingold(g1)

#lien A-D présent dans g1 et g2 donc g1-g2 et g2-g1 entrainent sa disparition
par(mfrow = c(2,2))
plot(g1, main = "g1", layout = lay, vertex.color = "red", vertex.size = 20, vertex.label.dist =3)
plot(g2, main = "g2", layout = lay, vertex.color = "red", vertex.size = 20, vertex.label.dist =3)
plot(diffg1, main = "graph.difference(g1, g2)", layout = lay, vertex.color = "red", vertex.size = 20, vertex.label.dist =3)
plot(diffg2, main = "graph.difference(g2, g1)", layout = lay, vertex.color = "red", vertex.size = 20, vertex.label.dist =3)
plot(diffg3, main = "graph.difference(g1, g3)", layout = lay, vertex.color = "red", vertex.size = 20, vertex.label.dist =3)

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






