#  created 09 Jan 2026
# Author Orchid ID: 0000-0001-7072-204X
# Author: Dr. Helmut C. Schmidt 
# helmut.schmidt@physics-beyond-standard-model.com
# Angaben gem. § 5 TMG

# Betreiber und Kontakt:
# Dr. Helmut Christian Schmidt
# Johann-Hackl-Ring 52
# D-85630
# Grasbrunn
# Germany
# Tel:+49 15256325349
# helmut.schmidt@physics-beyond-standard-model.com
# Sämtliche Inhalte auf dieser website (Texte, Titel, Bilder, Grafiken u. a.) unterliegen dem Schutz des Urheberrechts. 
# Sie dürfen außerhalb der Grenzen dieser Schutzgesetze nicht ohne die vorherige schriftliche Zustimmung vervielfältigt,
# verbreitet, veröffentlicht, verändert, Dritten zugänglich gemacht oder auf andere Weise genutzt werden.

# Verantwortlicher für journalistisch-redaktionelle Inhalte gem. § 55 II RstV:
# Dr. Helmut Christian Schmidt

import matplotlib.pyplot as plt  
import cmath
import matplotlib.colors as mcolors
pi = cmath.pi; i_T = 0; i_T1 = 0; E = 0; E = [0]*10;  g = [[0] * 10] * 10
Obj=[[0]*20]*30; N_E= [0]*20000; N_T= [0]*20000; Mark= [0]*20000; E_t=[0]*20000
Emax= [0]*1000; Emin= [0]*1000; i_Emax= [0]*30; i_Emin= [0]*30; Dmax = [[0]*10]*30; Dmin = [[0]*10]*30 

Obj=[["Name"   ,"m_e"               ,"E",  "-SD",  "+SD"                          ,"Charge","Spin","P.","name","pos"],
     ["e"      ,"1.00(1)"           ,"1.000","-0.005","0.005"                   ,"-1","1/2","","e",2],
     ["u"      ,"4.18(-0.51)(0.96)" ,"4.18","-0.51","0.96"                      ,"+2/3","","","u",3],
     ["d"      ,"9.14(-0.33)(0.94)" ,"9.14","-0.33","0.94"                      ,"-1/3","","","d",3],
     ["Muon"   ,"206.7682827(46)"   ,"206.7682827","-0.0000046","0.0000046"     ,"0","1","","muon",-6],
     ["s"      ,"182.8(-6.6)(16.8)" ,"182.8","-6.6","16.8"                      ,"-1/3","","","s",3],
     ["c"      ,"2485(-39)(39)"     ,"2485","-39","39"                          ,"+2/3","","","c",0],
     ["b"      ,"8186(14)"          ,"8186","-14","14"                          ,"-1/3","","","b",0],
     ["t"      ,"337710(570)"       ,"337710","-570","570"                      ,"+2/3","","","t",0],
     ["Higgs"  ,"244830(210)"       ,"244830","-210","210"                      ,"0","0","","Higgs",0],
     ["Pion 0" ,"264.1430(9)"       ,"264.1430","-0.0009","0.0009","0","0","-","$u\overline{d}-\overline{u}d$",3],        
     ["Pion +-","273.13243(35)"     ,"273.13243","-0.00035","0.00035","+-1","0","-",\
                                                                "$u\overline{u},\overline{d}d$",8], 
     ["Eta"    ,"1072.139(35)"      ,"1072.139","-0.035","0.035","0","0","-",\
                                                                "$u\overline{u}+\overline{d}d-2s\overline{s}$",1], 
     ["Eta`"   ,"1874.32(11)"       ,"1874.32","-0.11","0.11"   ,"0","0","-",\
                                                               "$u\overline{u}+\overline{d}d+s\overline{s}$",-22], 
     ["Rho +-" ,"1506(1)"           ,"1506","-1","1"            ,"-+1","1","-","$u\overline{u},\overline{d}d$",3], 
     ["Rho 0"  ,"1517.14(49)"       ,"1517.14","-0.49","0.49"   ,"0","1","-","$u\overline{u}-\overline{d}d$",8], 
     ["Omega"  ,"1531.62(25)"       ,"1531.62","-0.25","0.25"   ,"0","1","-","$u\overline{u}+\overline{d}d$",-8], 
     ["Phi"    ,"1995.035(31)"      ,"1995.035","-0.031","0.031","0","1","-","$s\overline{s}(most)$",-11],   
     ["K +-"   ,"966.102(21)"       ,"966.102","-0.021","0.021" ,"+-1","0","-","$u\overline{s},s\overline{u}$",2], 
     ["KL 0"   ,"973.800(26)"       ,"973.800","-0.026","0.026" ,"0","0","-","$d\overline{s},s\overline{d}$  ",7],  
     ["KS 0"   ,"973.800(26)"       ,"973.800","-0.026","0.026" ,"0","","","$d\overline{s},s\overline{d}$",12],   
     ["K* +-"  ,"1745.2(1)"         ,"1745.2","-0.1","0.1"      ,"+-1","","","$d\overline{s},s\overline{d}$",-12],
     ["K* 0"   ,"1752.6(1)"         ,"1752.6","-0.1","0.1"      ,"0","","","$d\overline{s},s\overline{d}$",-6],  
     ["Neutron","1838.68366200(74)" ,"1838.68366200", "-0.00000074", "0.00000074","0","1/2","1","udd",-5],  
     ["Proton" ,"1836.152673426(32)","1836.152673426","-0.000000032","0.000000032","1","1/2","1","uud",-10], 
     ["Tau"    ,"3477.23(23)"       ,"3477.23","-0.23","0.23"               ,"-1","1/2","","tau",0],    
     ["H"      ,"1837.47(-0.29)(0.20)","1837.47","-0.29","0.20"             ,"0","","","H",-5.5]]

def Energie(i4,i3,i2,i1,i0,i_1,C):
    g[2][4]=i4; g[2][3]=i3; g[2][2]=i2; g[1][1]=i1; g[1][0]=i0; g[1][-1]=i_1
    E[0]=0; E[1]=0; E[2]=0; E[3]=0; E[4]=0; E[5]=0; E[6]=0; E[7]=0
    
    E_C_pos = -pi + 2*pi**(-1) - pi**(-3) + 2*pi**(-5) - pi**(-7) + pi**(-9) - pi**(-12) - 2*pi**(-14) 
    E_C_neg = 2*pi - pi**(-1) + E_C_pos                                          
    E[0] = C*( (C>0) * E_C_pos - (C<0) * E_C_neg)                                    #Energy for the charge

    for l in range(4, 1, -1):                                                        #Gluonen r b g 
        E[2] += g[2][l]*(2*pi)**l
    for n in range(1, -2, -1):                                                       #e, u, d
        E[1] -= g[1][n]*(2*pi)**n
    
    for l in range(4, 1, -1):
        for n in range(1, -2, -1): 	     
            if g[2][l] != 0 and g[1][n] != 0:
                E[3] += (l+n<4)*(g[2][l]>0)* g[2][l] * g[1][n] * 2 * (2*pi)**(-l-n-1)#neutral, matter
                E[4] += (l+n<4)*(g[2][l]<0)* g[2][l] * g[1][n] * 2 * (2*pi)**(-l-n)  #neutral, antimatter
                E[5] -= (l+n>3) * g[2][l] * g[1][n] * 2 * (2*pi)**(-l-n-1)           #neutral, Gravitation
                E[6] += abs(g[2][l] * g[1][n]) * 2 * (2*pi)**(-8)                    #internal time
                g[2][l] = 0; g[1][n] = 0;  break
            if g[2][l] == 0 and g[1][n] == 0:
                E[7] -= (2*pi)**(-l-n-1)                                             #neutral, antimatter
                E[7] -= (2*pi)**(-l-n)                                               #neutral, antimatter   
                # -1/(2pi) >   Neutrino     \nu_{\mu} = 1/pi                          >  decay
                #              Antineutrino \nu_e = (2pi)**(-2(l+n)-1)/pi
                break                                                                   
    E[0] += E[1] + E[2] + E[3] + E[4] + E[5] + E[6] + E[7]
    return E[0]

for i5 in [0]:    # for speed v
    for i4 in [-1,-1/2,0,1/2,1]: #[-2,-1,-1/2,0,1/2,1,2]:  Select range  for more particle   for tau  [0]: #
        for i3 in [-1,-1/2,0,1/2,1]: #[-2,-1,-1/2,0,1/2,1,2]: # for Phi 
            for i2 in [-2,-3/2,-1,-1/2,0,1/2,1,3/2,2]:  
                for i1 in [-3,-2,-3/2,-1,-1/2,0,1/2,1,3/2,2,3]: #[-2,-3/2,-1,-1/2,0,1/2,1,3/2,2]:    
                    print("i4",i4,"i3",i3,"i2",i2,"i1",i1,"i_T1", i_T1)
                    for i0 in [-3,-2,-3/2,-1,-1/2,0,1/2,1,3/2,2,3]:
                        for i_1 in [-3,-2,-3/2,-1,-1/2,0,1/2,1,3/2,2,3]:  
                            for C in [-1,-1/2,0,1/2,1]:  
                                Energie(i4,i3,i2,i1,i0,i_1,C)
                                if E[0] > 2000: continue   # for plot  E < 2000
                                # Select any energy range. for H, Proton und Neutron: E[0] < 1836 or E[0] > 1839  
                                # if E[0] > 300: continue  for u, d, s, pion, muon
                                if  i4+i3+i2<0: continue # only E > 0
                                i_T += 1 
                                for j in [1,2,3,4,5,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26]: 
                                    min_ = float(Obj[j][3]);  max_ = float(Obj[j][4])
                                    if E[0]>0 and (E[0]-float(Obj[j][2])<=1.2*max_) and \
                                                                            (E[0]-float(Obj[j][2])>=1.2*min_):
                                        i_T1 += 1 
                                        N_E[i_T1]= E[0] 
                                        N_T[i_T1]= i_T  
                                        Mark[i_T1]= j 
                                        if Emax[j] < E[0]:
                                            Emax[j]= E[0]; i_Emax[j] = i_T; Dmax[j] = [i4,i3,i2,i1,i0,i_1,C]  
                                        if Emin[j] > E[0] or Emin[j] == 0:
                                            Emin[j]= E[0]; i_Emin[j] = i_T; Dmin[j] = [i4,i3,i2,i1,i0,i_1,C]  
                                if E[0] > 0:
                                    plt.plot([i_T,i_T+1], [E[0], E[0]], color="#C0BCBC")
print("..............   wait for Plot, several minutes   .................")
print("")
print("possible ET: ", i_T , "real ET: ",  i_T1)
F=["#FFFFFF","#000000","#F60000","#05FB4F","#000000","#CFCF00","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF",
   "#F700D2","#047619", "#00F73E","#7BB91F","#A9BF06","#789E20","#047619",
   "#E8EF88","#CF00B7","#EC61A9","#EF88BE","#88FF51","#FA9805",
   "#4200F6","#850827","#000000","#495999"]

for j in [1,2,3,4,5,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26]:
    m_min = float(Obj[j][2]) + float(Obj[j][3]); m_max = float(Obj[j][2]) + float(Obj[j][4])
    p = Obj[j][2].find('.') + len(str(int(m_min))) + 1
    E_mitte = (Emax[j]+Emin[j])/2; i_deltaE=abs(i_Emax[j]-i_Emin[j]+1); i_E = (i_Emax[j]+i_Emin[j])/2
    print("")
    print('{0:10}{1:16}{2:16}{3:8}{4:6}{5:5}{6:5}{7:5}{8:5}{9:5}{10:5}{3:16}'.format(Obj[j][0],"   meassure: m",\
                    "   theory: E"," total ","  i4","i3","i2","i1","i0","i-1","C",Obj[j][8]))
    print('{0:10}{1:16}{2:16}{3:8}{4:6}{5:5}{6:5}{7:5}{8:5}{9:5}{10:5}'.format("   max",round(m_max,p),\
                    round(Emax[j],p),i_Emax[j],Dmax[j][0],Dmax[j][1],Dmax[j][2],Dmax[j][3],Dmax[j][4],Dmax[j][5],Dmax[j][6]))    
    print('{0:10}{1:16}{2:16}{3:8}'.format("  mean",float(Obj[j][2]),round(E_mitte,8), i_deltaE))
    print('{0:10}{1:16}{2:16}{3:8}{4:6}{5:5}{6:5}{7:5}{8:5}{9:5}{10:5}'.format("   min",round(m_min,p),\
                    round(Emin[j],p),i_Emin[j],Dmin[j][0],Dmin[j][1],Dmin[j][2],Dmin[j][3],Dmin[j][4],Dmin[j][5],Dmin[j][6]))

    plt.text(i_E + 10000 * float(Obj[j][9]), E_mitte-20, Obj[j][8], color=F[j], fontsize=12)

for k in range(i_T1,0,-1):
    plt.plot(N_T[k], N_E[k], color=F[Mark[k]], marker='.', markerfacecolor=F[Mark[k]])  

   #plt.plot([1,i_T],[1836.1526734+1,1836.1526734+1],'k',linewidth=1);plt.text(1,1837,'m_Proton+1',\
   #                                                                                fontsize=12,color='blue')

   #  for 0 > E 2000
x_a= i_Emin[1] - 1000 ; x_e= i_T /2; plt.ylabel('Energy in $m_e$');  plt.xlim(i_Emin[1]-2000, i_T+30000); 
i2 = float(2*pi)**2; i3 = float(2*pi)**3; i4 = float(2*pi)**4; i5 = 1/2*(i4+i3+i2); i6 = i4+i3+i2;  
plt.plot([x_a,400000],[i4,i4],'k',linewidth=1);plt.text(x_a,i4+15,'$(2\pi)^4$', fontsize=12, color='blue')
plt.plot([x_a,400000],[i3,i3],'k',linewidth=1);plt.text(x_a,i3+15,'$(2\pi)^3$', fontsize=12, color='blue')
plt.plot([x_a,400000],[i2,i2],'k',linewidth=1);plt.text(x_a,i2+15,'$(2\pi)^2$', fontsize=12, color='blue')
plt.plot([x_a,400000],[i5,i5],'k',linewidth=1);plt.text(x_a,i5+15,'$1/2((2\pi)^4+(2\pi)^3+(2\pi)^2)$',\
                                                                                fontsize=12,color='blue')
plt.plot([x_a,400000],[i6,i6],'k',linewidth=1);plt.text(x_a,i6+15,'$(2\pi)^4+(2\pi)^3+(2\pi)^2$',\
                                                                                fontsize=12,color='blue')

x_a = i_T *3/4; dx = i_T/11; i = -50 
for j in [1,2,3,4,5,10,11,18,21,23,24,26]:
    particle = str(Obj[j][8]) 
    if j == 1  or j == 4  or j == 28:  particle=""          #  Obj[j][0].find(particle) == 0:
    plt.text(x_a, i, Obj[j][0]); plt.text(x_a+dx,i,particle); plt.text(x_a+2*dx, i, abs(i_Emax[j]-i_Emin[j]+1))
    i+= 80

plt.savefig('8.1.26.svg', format='svg', dpi=200)
plt.show()