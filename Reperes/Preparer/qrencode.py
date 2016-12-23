#!/usr/bin/python3

from PIL import Image
from sys import exit,stderr
from qrcorps import *
from qrcodestandard import *
from argparse import *

class qrencode:
  # Gestion des entrées (standard et options)

  def __init__(self):
    [self.arguments,self.nivcor,self.masque,self.message,self.mode,self.longclair,
      self.version,self.longtoutclair,self.clair,self.longueur,
      self.tout,self.entrelac,self.dim,self.tabmat,self.gris,self.bontab]=[None]*16


  def options(self):
    parser=ArgumentParser(description="Crée un code QR.")
    parser.add_argument("-m",choices=["0","1","2","3","4","5","6","7"],help="Force le choix du masque.")
    parser.add_argument("-i",required=True,metavar="fichier",help="Fichier texte d’entrée.")
    parser.add_argument("-c",metavar="image",help="Imagette placée au milieu du code.")
    parser.add_argument("-t",type=int,metavar="N",default=4,help="Taille d’un module.")
    parser.add_argument("-o",required=True,metavar="image",help="Image de sortie.")
    parser.add_argument("-n",choices=["L","M","Q","H"],default="L",help="Niveau de correction.")
    parser.add_argument("-a",choices=["0","1"],default="0",help="Afficher l’image ou non.")
    self.arguments=parser.parse_args()

  def entrees(self):
    if self.arguments is None:
      self.options()

    if self.arguments.t<=0:
      print("Il faut choisir une taille de module positive.")
      exit(1)
    self.nivcor=self.arguments.n
    if self.arguments.m is None:
      self.masque=None
    else:
      self.masque=int(self.arguments.m)
    self.message=""
    try:
      with open(self.arguments.i,"r") as f:
        for l in f:
          self.message=self.message+l
#      message=message[:-1]
    except IOError:
      print("Fichier %s inaccessible."%self.arguments.i)
      exit(1)
    #print(len(message))

# Détermination du mode et de la longueur des données

  def carac(self):
    if None in [self.nivcor,self.masque,self.message]:
      self.entrees()

    caracteres=set(self.message)

    if self.message.isdecimal():
      self.mode=0
      q,r=divmod(len(self.message),3)
      self.longclair=10*q
      if r==1:
        self.longclair+=4
      elif r==2:
        self.longclair+=7
    elif caracteres.issubset(set(alphanum)): # dans qrcodeoutils
      self.mode=1
      q,r=divmod(len(self.message),2)
      self.longclair=11*q+r*6
    else:
      self.mode=2
      self.message=self.message.encode("utf-8")
      self.longclair=len(self.message)*8

# Détermination de la version

  def detversion(self):
    if None in [self.mode,self.longclair]:
      self.carac()
    for self.version in tableau:
      self.longtoutclair=8*sum(blocs(tableau[self.version][self.nivcor])[1::2])
#      self.longtoutclair=8*sum(eval(tableau[self.version][self.nivcor].replace("×","*").replace("),",")+"))[1::2])
      s=self.longtoutclair-4-longbin(self.version,self.mode)
      if s>self.longclair:
        break

    if s<self.longclair:
      print("Le message est trop long, veuillez le raccourcir ou choisir un niveau de correction moins fort.")
      exit(1)

# Construction du tableau de données

  def donnees(self):
    if None in [self.longtoutclair,self.version]:
      self.detversion()
    self.clair=[0]*(3-self.mode)+[1]+[0]*self.mode
    self.longueur=dec2bin(len(self.message),longbin(self.version,self.mode))
    self.clair+=self.longueur

    if self.mode==0:
      q,r=divmod(len(self.message),3)
      for i in range(q):
        self.clair+=dec2bin(int(self.message[i*3:i*3+3]),10)
      if r==1:
        self.clair+=dec2bin(int(self.message[-1]),4)
      elif r==2:
        self.clair+=dec2bin(int(self.message[-2:]),7)
    elif self.mode==1:
      q,r=divmod(len(self.message),2)
      for i in range(q):
        self.clair+=dec2bin(45*alphanum.index(self.message[2*i])+alphanum.index(self.message[2*i+1]),11)
      if r==1:
        self.clair+=dec2bin(alphanum.index(self.message[-1]),6)
    else:
      for b in self.message:
        self.clair+=dec2bin(b,8)

# Remplissage par la fin de message et par le bourrage si nécessaire

  def fin(self):
    if None in [self.clair,self.longueur]:
      self.donnees()
    bourrage=dec2bin(0xec11,16)
    self.clair+=[0,0,0,0]
#    self.clair+=[0]*(8-len(self.clair)%8) # Pourtant, certains disent qu’il faut garder ça
    while len(self.clair)<self.longtoutclair:
      self.clair+=bourrage
    self.clair=self.clair[:self.longtoutclair]
#    print(len(clair))

# Calcul du code de Reed-Solomon

  def reedsolomon(self):
    if None in [self.clair,self.longueur]:
      self.fin()
    posclair,posredondant=court2blocs(tableau[self.version][self.nivcor])

#    longtout=8*sum(eval(tableau[version][nivcor].replace("×","*").replace("),",")+"))[::2])
#    longredondant=longtout-longclair
    blocsclair=[]
#    print(posclair,posredondant)
    for i in range(len(posclair)-1):
      blocsclair.append(self.clair[8*posclair[i]:8*posclair[i+1]])
    blocsredondant=[]
    for i in range(len(blocsclair)):
      bloc=blocsclair[i]
      longredondant=8*(posredondant[i+1]-posredondant[i])
      polyclair=message2poly(bloc+[0]*longredondant)
      modulo=Polynome.construction([1])
      for i in range(longredondant//8):
        modulo*=Polynome.construction([1,F256.exp(i)])
      blocsredondant.append(poly2message(polyclair%modulo))

#    print(len(blocsclair),len(blocsredondant))
#    print(blocsclair)
#    print(blocsredondant)
#    tout=clair+redondant
#    print(len(tout),version)
#    print(tout)
    self.tout=[]
    for bloc in blocsclair:
      self.tout+=bloc
#    print(len(tout))
    for bloc in blocsredondant:
      self.tout+=bloc
#    print(len(tout))

# Entrelacement du message binaire final

  def entrelacement(self):
    if self.tout is None:
      self.reedsolomon()

#    entrelac=blocs(tableau[version][nivcor])
#    for i in range(len(entrelac)//2):
#      entrelac[2*i]-=entrelac[2*i+1]

    posclair,posredondant=court2long(tableau[self.version][self.nivcor])
    pos=posclair+posredondant
    self.entrelac=[0]*len(self.tout)
    j=0
#    print(posclair)
    for i in range(len(self.entrelac)//8):
      if pos[i]!=-1:
        self.entrelac[8*j:8+8*j]=self.tout[8*pos[i]:8+8*pos[i]]
        j+=1
#    print(entrelac)
#    posclair,posredondant=court2blocs(tableau[version][nivcor])
#    clairblocs=[[] for _ in court2blocs(tableau[version][nivcor])[0][1:]] # court2blocs est dans qrcodestandard
#    redondantblocs=[[] for _ in court2blocs(tableau[version][nivcor])[0][1:]] # court2blocs vraiment utile ?
#    entrelac=[None for _ in range(8*max(posredondant))]
#    print(posclair,posredondant)
#    print(court2long(tableau[version][nivcor]))

#    for i in range(len(posclair)):
#      if posclair[i]!=-1:
#        entrelac[i*8:i*8+8]=clair[8*posclair[i]:8*(posclair[i]+1)]
#    print(entrelac)
#    for i in range(len(posredondant)):
#      if posredondant[i]!=-1:
#        entrelac[(max(posclair)+i+1)*8:(max(posclair)+i+1)*8+8]=redondantblocs[i%len(redondantblocs)][i%len(redondantblocs[i%len(redondantblocs)])]

# Création du tableau sans la zone de format

  def matrice(self):
    if self.entrelac is None:
      self.entrelacement()

    self.dim=17+4*self.version
    self.tabmat=[[0 for _ in range(self.dim)] for _ in range(self.dim)]
    for i in range(self.dim): # Échelles
      self.tabmat[6][i]=1-i%2
      self.tabmat[i][6]=1-i%2
    for i in range(7): # Cibles
      self.tabmat[i][:7]=cible[i]
      self.tabmat[-i-1][:7]=cible[i]
      self.tabmat[i][-7:]=cible[i]
    self.tabmat[-8][8]=1 # Le module toujours noir
    liste=minicibles[self.version] # Les minicibles
    for ci in liste:
      for cj in liste:
        if (ci,cj) not in {(6,6),(6,max(liste)),(max(liste),6)}:
          for i in range(ci-2,ci+3):
            self.tabmat[i][cj-2:cj+3]=minicible[i+2-ci]

# Remplissage du tableau

  def remplissage(self):
    if None in [self.dim,self.tabmat]:
      self.matrice()

    self.gris=griser(self.dim,self.version)
    i,j=self.dim-1,self.dim-1
    dire=-1
    k=0

    while k<len(self.entrelac):
      if self.gris[i][j]:
        self.tabmat[i][j]=self.entrelac[k]
        k+=1
      i,j,dire=suivant(i,j,dire,self.dim) # suivant est dans qrcodestandard

# Calcul du code de contrôle de la version

  def codecontrolev7(self):
    if None in [self.dim,self.tabmat]:
      self.remplissage()

    if self.version>=7:
      pol=[int(i) for i in "1111100100101"]
      liste=dec2bin(self.version,6)
      liste+=[0]*12
      while liste[0]==0:
        del(liste[0])
      while len(pol)<len(liste):
        pol.append(0)
      while len(liste)>12:
        liste=[i^j for (i,j) in zip(liste[:len(pol)],pol)]+liste[len(pol):]
        while liste[0]==0:
          del(liste[0])
      liste=[0]*(12-len(liste))+liste
      liste=dec2bin(self.version,6)+liste
      liste.reverse()

# Et placement de ce code

      for i in range(3):
        self.tabmat[-11+i][:6]=liste[i::3]
      for i in range(6):
        self.tabmat[i][-11:-8]=liste[3*i:3*i+3]

# Calcul du meilleur masque et placement des zones de formats

  def choixmasque(self):
    if None in [self.dim,self.tabmat]:
      self.codecontrolev7()

    cor=correctionniveau[self.nivcor]
    dmin=None

    for m in range(8):
      if self.masque is not None and self.masque!=m:
        continue
      tab=[list(l) for l in self.tabmat]
      ma=tuple(dec2bin(m,3))
      forma=cor+ma
      forma=forma+tuple(dec2bin(resteformat(bin2dec(forma)),10))
      forma=[i^j for (i,j) in zip(forma,masquef)]
      tab[8][:6]=forma[:6]
      tab[8][7:9]=forma[6:8]
      tab[7][8]=forma[8]
      for i in range(6):
        tab[i][8]=forma[-1-i]
      tab[8][-8:]=forma[-8:]
      for i in range(7):
        tab[-i-1][8]=forma[i]
      masquet=masques[ma]
      for i in range(self.dim):
        for j in range(self.dim):
          if self.gris[i][j]:
            tab[i][j]^=masquet(i,j)
      mal=malus(tab)
      if dmin is None or mal<dmin:
        dmin=mal
        self.bontab=tab
#      if mal<dmin:
#        self.bontab=tab
#        dmin=mal

# Création de l’image

  def creation(self):
    if self.bontab is None:
      self.choixmasque()

    self.taille=(8+self.dim)*self.arguments.t
    im=Image.new("RGB",(self.taille,self.taille),"white")
    blanc=(255,255,255)
    self.image=[blanc]*self.taille*4*int(self.arguments.t)
    for i in range(self.dim):
      for _ in range(self.arguments.t):
        self.image.extend([blanc]*4*self.arguments.t)
        for j in range(self.dim):
          couleur=255-255*self.bontab[i][j]
          self.image.extend([(couleur,couleur,couleur)]*self.arguments.t)
        self.image.extend([blanc]*4*self.arguments.t)
    self.image.extend([blanc]*self.taille*4*self.arguments.t)
    im.putdata(self.image)
    try:
      logo=Image.open(self.arguments.c)
      m=max(logo.width,logo.height)
      logo=logo.resize((self.taille//4*logo.width//m,self.taille//4*logo.height//m),Image.ANTIALIAS)
      im.paste(logo,(self.taille//2-logo.width//2,self.taille//2-logo.height//2),logo)
    except AttributeError:
      pass
    except IOError:
      print("Fichier %s inaccessible."%self.arguments.c)
      exit(1)

    if self.arguments.a=="1":
      im.show()
    im.save(self.arguments.o)

if __name__=="__main__":
  code=qrencode()
  code.creation()
