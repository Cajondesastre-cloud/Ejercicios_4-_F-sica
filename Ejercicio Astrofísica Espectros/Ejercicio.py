import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.io import fits
import scipy
from scipy.ndimage import filters

# lim se usa para dar un estimación de la distancia entre dos puntos contiguos o cercanos
# cuando se tiene un pico. De manera que será mayor cuanto más intenso sea el mismo y mayor
# flujo haya en el espectro.

# pw viene de peak wide, y de una idea del ancho promedio de los picos, de manera que se desprecian
# picos de amplitud menor. Esto es importante, puesto que dependiendo del espectro, la amplitud de 
# los mismos variará.

def iden(star, lines,labels, lim, pw):      # Para absorción.
    ident = []
    index_list = []
    for j in range(len(star[:,1])-pw):
        dif = star[j+pw,1] - star[j,1]
        if dif<-lim:
            for k in range(j,j+5):
                dif2 = star[k +pw,1] - star[k,1]
                if  dif2>lim:
                    index_list.append(k+int(pw/2))
                    break  
    picos = pd.unique(star[index_list[:],0])
    for i in range(len(picos)):
        for l in range(len(lines[:])):
            dif3 = picos[i] - lines[l]
            if abs(dif3)< 1.5:
                ident.append((lines[l], labels[l]))
    return ident

def iden_em(star, lines,labels, lim, pw):     # Para emisión.
    ident = []
    index_list = []
    for j in range(len(star[:,1])-pw):
        dif = star[j+pw,1] - star[j,1]
        if dif >lim:
            for k in range(j,j+5):
                dif2 = star[k +pw,1] - star[k,1]
                if  dif2<-lim:
                    index_list.append(k+int(pw/2))
                    break  
    picos = pd.unique(star[index_list[:],0])
    for i in range(len(picos)):
        for l in range(len(lines[:])):
            dif3 = picos[i] - lines[l]
            if abs(dif3)< 2:
                ident.append((lines[l], labels[l]))
    return ident

def color_def(label):       # Color dependiendo de la etiqueta.
    c = "white"
    if label == "CII":
        c = "red"
    elif label == "CaII":
        c = "brown"
    elif label == "TiI":
        c = "yellow"
    elif label == "TiII":
        c = "orange"
    elif label == "SiIII":
        c = "purple"
    elif label == "MgII":
        c = "pink"
    elif label == "H$\\epsilon$":
        c = "green"
    elif label == "H$\\delta$":
        c = "limegreen"
    elif label == "H$\\gamma$":
        c = "lime"
    elif label == "H$\\beta$":
        c = "springgreen"
    elif label == "HeI":
        c = "blue"
    elif label == "HeII":
        c = "royalblue"
    elif label == "VI":
        c = "black"
    elif label == "NII":
        c = "aquamarine"
    elif label == "NIII":
        c = "turquoise"
    elif label == "NIV":
        c = "lightseagreen"
    elif label == "NV":
        c = "teal"
    elif label == "SiIV":
        c = "fuchsia"
    elif label =="MnI":
        c ="olive"
    elif label == "SiII":
        c = "mediumpurple"
    elif label == "TiO":
        c = "orangered"
    elif label =="FeI":
        c = "darkslategrey"
    
    return c

# Identificaciones

star_1 = np.loadtxt("starprob1.dat")
#A9V = np.loadtxt("s0921")
A6IV_V = np.loadtxt("s0920")

star_2 = np.loadtxt("starprob2.dat")
B3IV = np.loadtxt("s0638")
B3III = np.loadtxt("s0817")
B1 = np.loadtxt("s0737")
#B2_5IV = np.loadtxt("s0715")
B3Ia = np.loadtxt("s0778")

star_3 = np.loadtxt("starprob3.dat")
w=np.genfromtxt('starprob3.dat',usecols=(0))
f=np.genfromtxt('starprob3.dat',usecols=(1))
smoothed=scipy.ndimage.filters.uniform_filter(f,size=15)
M2Iabs = np.loadtxt("s0186")

star_4 = np.loadtxt("starprob4.dat")
# Abrimos el espectro de la estrella que queramos en formato fits
file=fits.open("_wr6_20111021_103.fits")
# La siguiente linea recupera los datos de las longitudes de onda
wave = file[0].header['CRVAL1'] + np.arange(len(file[0].data))*file[0].header['CDELT1']
#  La siguiente linea recupera los datos del flujo
flux=file[0].data

# Identificación de líneas.

lineas = np.loadtxt("lineas.txt", usecols=1)
labels = np.loadtxt("lineas.txt", dtype=str, usecols=0)
lin_ley = np.loadtxt("Leyenda.txt", usecols=1)
lab_ley = np.loadtxt("Leyenda.txt", dtype=str, usecols=0)

ident_1 = iden(star_1, lineas, labels, 70, 25)
ident_2 = iden(star_2, lineas, labels, 100, 7)
ident_4_abs = iden(star_4, lineas, labels, 10, 2)
ident_4_em = iden_em(star_4, lineas, labels, 12, 50)

#Plots

fig, axs = plt.subplots(2, 2)
plt.suptitle("Clasificación espectral de estrellas problema")

axs[0,0].set_xlabel('$\lambda [\mathrm{\AA}$]')
axs[0,0].set_ylabel("F($\lambda$)")
axs[0,0].set_xlim([star_1[0, 0],star_1[-1, 0]])
axs[0,0].plot(star_1[:, 0], star_1[:,1], color="black")
for i in range(len(ident_1)):
    axs[0,0].vlines(ident_1[i][0], 0, 5000, color=color_def(ident_1[i][1]), linestyle="dashed", alpha = 0.7, label=ident_1[i][1])
    axs[0,0].annotate(ident_1[i][1], (ident_1[i][0],1000 +500*(-1)**i), textcoords="offset points",xytext=(10,10), 
                 ha='center', color=color_def(ident_1[i][1]))
axs[0,0].plot(A6IV_V[:, 0], A6IV_V[:,1]/max(A6IV_V[:,1])*max(star_1[:,1]), color="purple", alpha=0.6)
axs[0,0].vlines(4128, 0, 5000, color=color_def("SiII"), linestyle="dashed", alpha = 0.7, label="SiII")
axs[0,0].annotate("SiII", (4128,1000), textcoords="offset points",xytext=(10,10), 
                 ha='center', color=color_def("SiII"))

axs[0,1].set_xlabel('$\lambda [\mathrm{\AA}$]')
axs[0,1].set_ylabel("F($\lambda$)")
axs[0,1].set_xlim([star_2[0, 0],star_2[-1, 0]])
axs[0,1].plot(star_2[:, 0], star_2[:, 1],color="black")
for i in range(len(ident_2)):
    axs[0,1].vlines(ident_2[i][0], 0, 5000, color=color_def(ident_2[i][1]), linestyle="dashed", alpha = 0.7, label=ident_2[i][1])
    axs[0,1].annotate(ident_2[i][1], (ident_2[i][0],1500+500*(-1)**i), textcoords="offset points",xytext=(10,10), 
                 ha='center', color=color_def(ident_2[i][1]))
axs[0,1].plot(B1[:, 0], B1[:, 1]/max(B1[:,1])*max(star_2[:,1])*0.8, color="purple", alpha=0.6)
axs[0,1].vlines(4089, 0, 5000, color=color_def("SiIV"), linestyle="dashed", alpha = 0.7, label="SiIV")
axs[0,1].annotate("SiIV", (4089,1000), textcoords="offset points",xytext=(10,10), 
                 ha='center', color=color_def("SiIV"))


axs[1,0].set_xlabel('$\lambda [\AA}$]')
axs[1,0].set_ylabel("F($\lambda$)")
axs[1,0].set_xlim([star_3[0, 0],star_3[-1, 0]])
axs[1,0].plot(star_3[:, 0], smoothed[:],color="black")
axs[1,0].plot(M2Iabs[:, 0], M2Iabs[:, 1]/max(M2Iabs[:,1])*max(star_3[:,1])*1.1, color="purple", alpha=0.6)
axs[1,0].vlines(4584, 0, 60000, color=color_def("TiO"), linestyle="dashed", alpha = 0.7, label="TiO")
axs[1,0].annotate("TiO", (4584,50000), textcoords="offset points",xytext=(-10,10), 
                 ha='center', color=color_def("TiO"))
axs[1,0].vlines(4623, 0, 60000, color=color_def("TiO"), linestyle="dashed", alpha = 0.7, label="TiO")
axs[1,0].vlines(4759, 0, 60000, color=color_def("TiO"), linestyle="dashed", alpha = 0.7, label="TiO")
axs[1,0].vlines(4804, 0, 60000, color=color_def("TiO"), linestyle="dashed", alpha = 0.7, label="TiO")
axs[1,0].vlines(4954, 0, 60000, color=color_def("TiO"), linestyle="dashed", alpha = 0.7, label="TiO")


axs[1,1].set_xlabel('$\lambda [\mathrm{\AA}$]')
axs[1,1].set_ylabel("F($\lambda$)")
axs[1,1].set_xlim([star_4[0, 0],star_4[-1, 0]])
axs[1,1].plot(star_4[:, 0], star_4[:, 1], color="black")
axs[1,1].plot(wave[:], flux[:]/max(flux[:])*max(star_4[:,1]), color="purple", alpha=0.6)
for i in range(len(ident_4_abs)):
    axs[1,1].vlines(ident_4_abs[i][0], 0, 6000, color=color_def(ident_4_abs[i][1]), linestyle="dashed", alpha = 0.7, label=ident_4_abs[i][1])
    axs[1,1].annotate(ident_4_abs[i][1], (ident_4_abs[i][0],1500+500*(-1)**i), textcoords="offset points",xytext=(10,10), 
                 ha='center', color=color_def(ident_4_abs[i][1]))
for i in range(len(ident_4_em)):
    axs[1,1].vlines(ident_4_em[i][0], 0, 6000, color=color_def(ident_4_em[i][1]), linestyle="dashed", alpha = 0.7, label=ident_4_em[i][1])
    axs[1,1].annotate(ident_4_em[i][1], (ident_4_em[i][0],1500+500*(-1)**i), textcoords="offset points",xytext=(10,10), 
                 ha='center', color=color_def(ident_4_em[i][1]))
axs[1,1].vlines(4606.33, 0, 6000, color=color_def("NIV"), linestyle="dashed", alpha = 0.7, label="NIV")
axs[1,1].annotate("NIV", (4606.33,1000), textcoords="offset points",xytext=(10,10), 
                 ha='center', color=color_def("NIV"))
axs[1,1].vlines(4944.56, 0, 6000, color=color_def("NV"), linestyle="dashed", alpha = 0.7, label="NV")
axs[1,1].annotate("NV", (4944.56,2000), textcoords="offset points",xytext=(10,10), 
                 ha='center', color=color_def("NV"))

plt.tight_layout()
plt.savefig('Espectros.pdf')
plt.show()


Star_1 = plt.figure()
plt.xlabel('$\lambda [\mathrm{\AA}$]')
plt.ylabel("F($\lambda$)")
plt.xlim([star_1[0, 0],star_1[-1, 0]])
plt.plot(star_1[:, 0], star_1[:,1], color="black")
for i in range(len(ident_1)):
    plt.vlines(ident_1[i][0], 0, 5000, color=color_def(ident_1[i][1]), linestyle="dashed", alpha = 0.7, label=ident_1[i][1])
    plt.annotate(ident_1[i][1], (ident_1[i][0],1000 +500*(-1)**i), textcoords="offset points",xytext=(10,10), 
                 ha='center', color=color_def(ident_1[i][1]))
plt.plot(A6IV_V[:, 0], A6IV_V[:,1]/max(A6IV_V[:,1])*max(star_1[:,1]), color="purple", alpha=0.6)
plt.vlines(4128, 0, 5000, color=color_def("SiII"), linestyle="dashed", alpha = 0.7, label="SiII")
plt.annotate("SiII", (4128,1000), textcoords="offset points",xytext=(10,10), 
                 ha='center', color=color_def("SiII"))
plt.show()


Star_2 = plt.figure()
plt.xlabel('$\lambda [\mathrm{\AA}$]')
plt.ylabel("F($\lambda$)")
plt.xlim([star_2[0, 0],star_2[-1, 0]])
plt.plot(star_2[:, 0], star_2[:, 1],color="black")
for i in range(len(ident_2)):
    plt.vlines(ident_2[i][0], 0, 5000, color=color_def(ident_2[i][1]), linestyle="dashed", alpha = 0.7, label=ident_2[i][1])
    plt.annotate(ident_2[i][1], (ident_2[i][0],1500+500*(-1)**i), textcoords="offset points",xytext=(10,10), 
                 ha='center', color=color_def(ident_2[i][1]))
plt.plot(B1[:, 0], B1[:, 1]/max(B1[:,1])*max(star_2[:,1])*0.8, color="purple", alpha=0.6)
plt.vlines(4089, 0, 5000, color=color_def("SiIV"), linestyle="dashed", alpha = 0.7, label="SiIV")
plt.annotate("SiIV", (4089,1000), textcoords="offset points",xytext=(10,10), 
                 ha='center', color=color_def("SiIV"))
plt.show()


Star_3 = plt.figure()
plt.xlabel('$\lambda [\AA}$]')
plt.ylabel("F($\lambda$)")
plt.xlim([star_3[0, 0],star_3[-1, 0]])
plt.plot(star_3[:, 0], smoothed[:],color="black")
plt.plot(M2Iabs[:, 0], M2Iabs[:, 1]/max(M2Iabs[:,1])*max(star_3[:,1])*1.1, color="purple", alpha=0.6)
plt.vlines(4584, 0, 60000, color=color_def("TiO"), linestyle="dashed", alpha = 0.7, label="TiO")
plt.annotate("TiO", (4584,50000), textcoords="offset points",xytext=(-10,10), 
                 ha='center', color=color_def("TiO"))
plt.vlines(4623, 0, 60000, color=color_def("TiO"), linestyle="dashed", alpha = 0.7, label="TiO")
plt.vlines(4759, 0, 60000, color=color_def("TiO"), linestyle="dashed", alpha = 0.7, label="TiO")
plt.vlines(4804, 0, 60000, color=color_def("TiO"), linestyle="dashed", alpha = 0.7, label="TiO")
plt.vlines(4954, 0, 60000, color=color_def("TiO"), linestyle="dashed", alpha = 0.7, label="TiO")
plt.show()

Star_4 = plt.figure()
plt.xlabel('$\lambda [\mathrm{\AA}$]')
plt.ylabel("F($\lambda$)")
plt.xlim([star_4[0, 0],star_4[-1, 0]])
plt.plot(star_4[:, 0], star_4[:, 1], color="black")
plt.plot(wave[:], flux[:]/max(flux[:])*max(star_4[:,1]), color="purple", alpha=0.6)
for i in range(len(ident_4_abs)):
    plt.vlines(ident_4_abs[i][0], 0, 6000, color=color_def(ident_4_abs[i][1]), linestyle="dashed", alpha = 0.7, label=ident_4_abs[i][1])
    plt.annotate(ident_4_abs[i][1], (ident_4_abs[i][0],1500+500*(-1)**i), textcoords="offset points",xytext=(10,10), 
                 ha='center', color=color_def(ident_4_abs[i][1]))
for i in range(len(ident_4_em)):
    plt.vlines(ident_4_em[i][0], 0, 6000, color=color_def(ident_4_em[i][1]), linestyle="dashed", alpha = 0.7, label=ident_4_em[i][1])
    plt.annotate(ident_4_em[i][1], (ident_4_em[i][0],1500+500*(-1)**i), textcoords="offset points",xytext=(10,10), 
                 ha='center', color=color_def(ident_4_em[i][1]))
plt.vlines(4606.33, 0, 6000, color=color_def("NIV"), linestyle="dashed", alpha = 0.7, label="NIV")
plt.annotate("NIV", (4606.33,1000), textcoords="offset points",xytext=(10,10), 
                 ha='center', color=color_def("NIV"))
plt.vlines(4944.56, 0, 6000, color=color_def("NV"), linestyle="dashed", alpha = 0.7, label="NV")
plt.annotate("NV", (4944.56,2000), textcoords="offset points",xytext=(10,10), 
                 ha='center', color=color_def("NV"))
plt.show()

fig2 = plt.figure()
plt.xlim([3800,5500])
for i in range(len(lin_ley)):
    plt.vlines(lin_ley[i], 0, 4000, colors = color_def(lab_ley[i]), label=lab_ley[i], linestyle="dashed")
plt.legend(loc="best")
plt.savefig("Leyenda.pdf")
