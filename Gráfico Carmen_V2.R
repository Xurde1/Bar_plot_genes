

#################################################################################################################################
################################################ Gr?fico de barras para Carmen ##################################################
#################################################################################################################################


#El script es para generar unos gr?ficos con los genes de cada uno de los enrchiments realizados previamente

#el input es un txt con los datos necesarios y se usa ggplot2 para la generaci?n de gr?ficos

library(ggplot2)

#Indentificar el directorio de trabajo//carpeta donde estan los ficheros .txt
#setwd("~/Sexo psoriasis/GSE63315_UV/Common/")

#Indicar el archivo txt de entrada
archivo<-"males DSigDB_table.txt"

#carga el archivo
input<-read.delim(archivo, header = TRUE, sep = "\t", dec = ".")

#ordenamos por FDR
input<-input[order(input$Adjusted.P.value), ]


#formatear el data frame

#n?mero de filas que queremos conservar; s?lo las que tiene FDR < 0.05
trues<-input$Adjusted.P.value < 0.05        
casos<-length(trues[trues ==TRUE])

#Si no hay  gene sets significativos se toman los 3 primeros
if (casos ==0){
  casos = 3
  
}else casos =casos



#Si hay m?s de 5 gene sets significativos se toman solo los 5 primeros
if (casos >5){
  casos = 5

}else casos =casos


#eliminanos las columnas no necesarias y nos quedamos con los casos que queremos

datos<-input[1:casos, c(1,4, 7, 8) ]
genes<-input[1:casos,c(1,9)]

#contamos el n?mero de genes que aparece en cada gene set
num_genes<-lengths(gregexpr(";", genes$Genes)) + 1

#creamos la tabla con los genes contados
tabla<-cbind(datos, num_genes)

#ID para la tabla
ID<-seq.int(1,casos, 1)
tabla<-cbind(ID, tabla)#creamos un ID para cada fila, esto es porque los Term son muy largos y no se pueden poner en el eje x

#ordenamos la tabla por el FDR
tabla<-tabla[order(tabla$Adjusted.P.value), ]


#Tomaos por separado las listas de genes de cada gene set
lista_genes<-input[1:casos, ]
vector_genes<-as.vector(lista_genes[,9])


################################ si tenemos m?s de 40 genes y queremos representar 40 de cada gene set  #########################

#por si tenemos m?s de 40 genes en el gene set; el en gr?fico solo entran 40 genes como m?ximo

#vector_genes2<-list()
#con esto se representan 40 genes de cada gene set
#for (i in 1:casos){
  
  #gene_pos<-unlist(gregexpr(pattern = ";", vector_genes[i]))
  #gene_pos[40]
  #vector_genes2[i]<-substr(vector_genes[i], 0, gene_pos[40]-1)
#}

#################################################################################################################################

#El m?ximo n?mero de genes que se puede representar es 40
maximo<-40

#En esta condicional si tiene mas genes que elmaximo se reduce de forma proporcional en numero de genes

if (max(tabla$num_genes)> maximo){
  #para representar un n?mero proporcional de genes de cada gene set
  cuantos_genes<-as.vector(tabla$num_genes)
  genes_proporcional<-list()
  porcentage_representado<-list()
  gene_pos<-list()
  
  
  for (i in 1:casos) {
    genes_proporcional[i]<-round((maximo * cuantos_genes[i])/max(cuantos_genes), 0) #numero de genes a representar
    porcentage_representado[i]<-round(((maximo * 100)/cuantos_genes[i]), 2) #obtenemos el porcentage de gebnes representados por gene set 
    
    gene_pos<-unlist(gregexpr(pattern = ";", vector_genes[i]))
    vector_genes[i]<-substr(vector_genes[i], 0, gene_pos[as.numeric(genes_proporcional[i])]-1)
  }
  tabla = cbind(tabla, as.numeric(porcentage_representado))
}else vector_genes = vector_genes





#como tenemos el eje x con un factor y ggplot los ordenano alfab?ticamente por defecto, especificamos el orden deseado
tabla$Term <-factor(tabla$Term, levels = tabla$Term[order(tabla$Adjusted.P.value)])


#listas para anotar en el gr?fico

new_list<-list() #lista para guardar los genes en un a columna
nueva_lista<-list()#lista para dibujar la lista de genes en el gr?fico
Anotaciones<-list()#lista para guardar los scores
Colores<-list("#000000", "#808080", "#A9A9A9", "#C0C0C0", "#D3D3D3") #esacala de grises(5)

#limite de altura que tendr? el eje y
limit<-round(max(tabla$num_genes)+(max(tabla$num_genes)*0.2), 0)




#Bucle para crear todas las anotaciones el gr?fico ;agrupando las anotaciones
for (i in 1:casos) {
  new_list[i]<- gsub(";", "\n", vector_genes[i])
  
  nueva_lista[[i]]<- annotate("text", x= i , y=0, colour=Colores[[i]], label= new_list[[i]], size=1.9, vjust= 0, fontface =2)
  
  Anotaciones[[i]]<-annotate("text", x=i+0.2, y=max(tabla$num_genes[[i]])+5, paste0("OR: ", round(tabla$Odds.Ratio[[i]], 3)), paste0("CS: ", round(tabla$Combined.Score[[i]], 1)),size= 3, fontface =2, hjust=0.05)
                                                                                       

}



#Gr?fico

title=gsub("_table.txt", "", paste0("Represented gene sets in ",archivo)) 

plot<-ggplot(tabla, aes(x=ID, y=num_genes))

plot+
  geom_point(color="transparent")+
  nueva_lista +
  Anotaciones+
  labs(title=title, x="Gene Sets", y="Number of genes")+#, caption = "Produced by Hospital la Princesa")+
  ylim(0, limit)+
  xlim(0.5, casos+0.5)+
  theme(axis.text.y   = element_text(size=14),
        axis.text.x   = element_text(size=14),
        axis.title.y  = element_text(size=14),
        axis.title.x  = element_text(size=14),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)
  )


#en esta tabla se resumen los datos del grafico  
names(tabla)[7]="percentage of genes in graphic"