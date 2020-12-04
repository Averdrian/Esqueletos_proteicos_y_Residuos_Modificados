import re

def name(lines):
    buscados = []
    for i in lines:
        if re.search('^TITLE', i):
            buscados.append(i)
    for i in range(len(buscados)):
        #elimina la palabra TITLE y tambien los espacios en blanco adicionales a la izquierda y derecha junto con la numeracion de linea

        buscados[i] = re.sub("^TITLE", "", buscados[i])
        buscados[i] = re.sub("$(A-Za-z)*", "", buscados[i])
        buscados[i] = buscados[i].strip()
        buscados[i] = re.sub("^([0-9])*", "" , buscados[i])
    nombre = "".join(buscados)
    return nombre

def residuosModificados(lines):
    grupos = []
    gruposFinal = []
    for i in range(len(lines)):
        if re.search("^MODRES", lines[i]):grupos.append(lines[i])
    #para cada uno le cambiamos el formato y si el residuo no esta aun en la lista se añade a ella
    for i in range(len(grupos)):
        grupos[i] = re.sub("^MODRES ", "", grupos[i]).strip()
        grupos[i] = re.sub("^[1-9A-Z]* ", "", grupos[i])
        aux = re.search("^[A-Z]*", grupos[i]).group()
        if not gruposFinal.__contains__(aux): gruposFinal.append(aux)
    return gruposFinal



def esqueletoProteico(lines):
    buscados = []
    for i in lines:
        if re.search("^(ATOM|HETATM)", i):
            buscados.append(i)
    for i in range(len(buscados)):
        buscados[i] = re.sub("^(ATOM|HETATM)", "", buscados[i])
        buscados[i] = buscados[i].strip()
        buscados[i] = re.sub("^([0-9])*", "", buscados[i])
        buscados[i] = re.sub("([A-Z]*)$", "", buscados[i])
        buscados[i] = buscados[i].strip()

    agrupaciones = []
    #agrupaciones contendra la letra identificadora de cada cadena
    for i in buscados:
        esqueletos = re.search(".* ([A-Z]) .*", i)
        i = re.sub(".* ([A-Z]) .*", "", i)
        if esqueletos.group(1) not in agrupaciones: agrupaciones.append(esqueletos.group(1))
    for i in agrupaciones:
        fichSal = open('5ujw-'+i+'-.txt', 'w')
        fichSal.write("".ljust(12, " ") + "N".ljust(24, " ") + "CA".ljust(24, " ") + "C".ljust(24, " ") + "O".ljust(24, " ") + "\n")
        groups = []
        #busca los que tengan los atomos requeridos y los guarda en una lista
        for j in range(len(buscados)):
            if re.search("^((C|CA|N|O) )", buscados[j]) \
                    and re.search(" {} ".format(i), buscados[j]) \
                    and not re.search(" HOH ", buscados[j]): groups.append(buscados[j])

        #cogera los grupos de cuatro en cuatro (uno por cada atomo) y los procesara
        while(len(groups) > 0):

            subgrupo = groups[0:4]
            groups = groups[4:]
            num = re.search(" ([0-9]+) ", subgrupo[1]).group(1)
            prot = re.search(".*([A-Z]{3})", subgrupo[1]).group(1)
            for y in range(len(subgrupo)):
                subgrupo[y] = re.sub(" {} ".format(num), "", subgrupo[y])
                subgrupo[y] = re.sub("{}".format(prot), "", subgrupo[y])
                subgrupo[y] = re.sub("^[A-Z]*", "", subgrupo[y]).strip()
                subgrupo[y] = re.sub(" [A-Z] "," ",subgrupo[y])
                subgrupo[y] = re.sub("1\.00.*", "", subgrupo[y])
                subgrupo[y] = re.sub("^[A-Z]","",subgrupo[y]).strip()
                subgrupo[y] = subgrupo[y].ljust(22, " ")
                s = subgrupo[y].split(" ")
                for z in s:
                    z.rjust(7, " ")
                subgrupo[y] = " ".join(s)
            #escribimos los datos requeridos en el fichero
            fichSal.write(num.ljust(5, " ") + "  " + prot + "  " + "  ".join(subgrupo) + "\n")

        #al terminar con cada fichero lo cerramos
        fichSal.close()



def noPeptidicas(lines, modres):
    l = []
    #en la lista l se guardaran las lineas validas (etiqueta HETATM que no esten en MODRES y no sean agua)
    for i in lines:
        if re.search("^HETATM", i):
            for j in modres:
                if not re.search(j, i) and not re.search("HOH", i): l.append(i)

    #para cada linea le cambiamos el formato al pedido en el ejercicio y lo añadimos a su correspondiente archivo
    for i in range(len(l)):
        #todas las llamadas a strip() son para suprimir espacios en blanco en los extremos de la cadena
        #quitamos la primera palabra que haya en cadena
        l[i] = re.sub("^[A-Z0-9]*", "", l[i]).strip()
        #quitamos la palabra que haya al final de la cadena
        l[i] = re.sub("([A-Z]*)$", "", l[i]).strip()
        #quitamos la otra palabra innecesaria tambien al inicio de la cadena
        molecula = re.sub("^[A-Z0-9]*", "", l[i]).strip()
        #nos quedamos solo con la que ahora este al inicio de la cadena y possteriormente borramos dicha palabra del elemento de la lista
        molecula = re.search("^[A-Z0-9]*", molecula).group()
        l[i] = re.sub(molecula, "", l[i])
        #la letra que identifica a que grupo la guardamos para poder usarla al abrir el archivo de salida
        #despues tambien la borramos, es la unica letra aislada (con espacios en blanco a ambos lados) que puede quedar en la cadena
        f = re.search(" [A-Z] ", l[i]).group().strip()
        l[i] = re.sub(f, "", l[i])
        #todas las lineas tienen dos numeros que no son requeridos, uno de ellos como se puede ver un '1.00' y otro un natural
        l[i] = re.sub("1\.00.*", "", l[i]).strip()
        #el natural que buscamos sabemos que puede tener de 1 a n caracteres (nunca 0), por eso usamos '+'
        l[i] = re.sub(" [0-9]+ ", "", l[i])
        #el archivo se abre en modo 'a' porque varias iteraciones van a modificar el mismo archivo
        fh = open("5ujw-" + f + "-" + molecula + ".txt", 'a')
        fh.write(l[i] + "\n")
        #al final siempre cerrar el archivo
        fh.close()




fich = open('5ujw.pdb', 'r')
lines = fich.readlines()
fich.close()
nombre = name(lines)
esqueletoProteico(lines)
resModificados = residuosModificados(lines)
noPeptidicas(lines, resModificados)
print(nombre)
print("Residuos modificados:", resModificados)


