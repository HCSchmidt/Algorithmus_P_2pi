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
import numpy as np
Op = 2; 
if Op==1: O0=1; O1=1; f_="Polynom 1113331 a 26.1.26.txt"; f_png = "Polynom 111333 1 a 26.1.26.png"  # fast ,
if Op==2: O0=2; O1=2; f_="Polynom 2223331 a 26.1.26.txt"; f_png = "Polynom 2223331 a 26.1.26.png" # slow , more particles
if Op==3: O0=0; O1=1; f_="Polynom 0113331 a 26.1.26.txt"; f_png = "Polynom 0113331 a 26.1.26.png" # fast, for u, d, s, pion, muon E<300
if Op==4: O0=1; O1=1; f_="Polynom 1113331 Atom a 26.1.26.txt"; f_png = 'Polynom 1113331 Atom a 26.1.26.png' # fast, for Proton H-Atom Neutron
if Op==5: O0=2; O1=2; f_="Polynom 2223331 Tau a 26.1.26.txt"; f_png = "Polynom 2223331 Tau a 26.1.26.png" # slow, more particles with c, Deuterium,tau

f=open(f_,'w', encoding="utf8"); f=open(f_,'a+', encoding="utf8")
pi = cmath.pi; i_T = 0; i_T1 = 0; E = 0; i1=0;i0=0;i_1=0; E = [0]*10; m=0; g = [[0]*10]*10; j=0; D_i_E = [0]*300;D_i_N=[0]*300
Obj=[[0]*20]*30; N_E= [0]*100000; N_T= [[0]*100]*100; E_t=[0]*100000; X=[0]*30
Emax = np.zeros((30, 513), dtype=float); Emin = np.zeros((30, 513), dtype=float)
i_Emax = np.zeros((30, 513), dtype=int); i_Emin = np.zeros((30, 513), dtype=int)
Dmax = np.zeros((30, 513, 7), dtype=int); Dmin = np.zeros((30, 513, 7), dtype=int)

Obj=[["Name"   ,"m_e"               ,"E",   "-SD",   "+SD"                      ,"Halbwertszeit T in sec" ,"Charge","Spin","P.","name","pos"],
     ["e"      ,"1.00000000000(31)" ,"1.00000000000","-0.005","0.000"           ," "                       ,"-1","1/2","","e",3],  
     ["u"      ,"4.18(-0.51)(0.96)" ,"4.18","-0.51","0.96"                      ," "                       ,"+2/3","","","u",6],
     ["d"      ,"9.14(-0.33)(0.94)" ,"9.14","-0.33","0.94"                      ," "                       ,"-1/3","","","d",9],
     ["s"      ,"182.8(-6.6)(16.8)" ,"182.8","-6.6","16.8"                      ," "                       ,"-1/3","","","s",3],
     ["Muon"   ,"206.7682827(46)"   ,"206.7682827","-0.0000046","0.0000046"     ,"2.1969811(22) e -6"     ,"0","1","","muon",4],
     ["Pion 0" ,"264.1430(9)"       ,"264.1430","-0.0009","0.0009"              ,"8.52(18) e -17"     ,"0","0","-","$u\overline{d}-\overline{u}d$",4],        
     ["Pion +-","273.13243(35)"     ,"273.13243","-0.00035","0.00035"           ,"2.6033(5) e -8"     ,"+-1","0","-","$u\overline{u},\overline{d}d$",8], 
     ["K +-"   ,"966.102(21)"       ,"966.102","-0.021","0.021"                 ,"1.2380(20) e -8"     ,"+-1","0","-","$u\overline{s},s\overline{u}$",-8], 
     ["KL 0"   ,"973.800(26)"       ,"973.800","-0.026","0.026"                 ,"5.116(21) e -8"      ,"0","0","-","$d\overline{s},s\overline{d}$  ",-16],  
     ["KS 0"   ,"973.800(26)"       ,"973.800","-0.026","0.026"                 ,"8.954(4) e -11"      ,"0","","","$d\overline{s},s\overline{d}$",-8],   
     ["Eta"    ,"1072.139(35)"      ,"1072.139","-0.035","0.035"                ,"5 e -19"             ,"0","0","-","$u\overline{u}+\overline{d}d-2s\overline{s}$",-15], 
     ["Rho +-" ,"1506(1)"           ,"1506","-1","1"                            ,"4 e -24"             ,"-+1","1","-","$u\overline{u},\overline{d}d$",3], 
     ["Rho 0"  ,"1517.14(49)"       ,"1517.14","-0.49","0.49"                   ,"4 e -24"             ,"0","1","-","$u\overline{u}-\overline{d}d$",8], 
     ["Omega"  ,"1531.62(25)"       ,"1531.62","-0.25","0.25"                   ,"7.75(7) e -23"       ,"0","1","-","$u\overline{u}+\overline{d}d$",-8], 
     ["K* +-"  ,"1745.2(1)"         ,"1745.2","-0.1","0.1"                      ,"1.3 e -23"           ,"+-1","","","$d\overline{s},s\overline{d}$",-16],
     ["K* 0"   ,"1752.6(1)"         ,"1752.6","-0.1","0.1"                      ,"1.3 e -23"           ,"0","","","$d\overline{s},s\overline{d}$",-6],  
     ["Proton" ,"1836.152673426(32)","1836.152673426","-0.000000032","0.000000032"," "                 ,"1","1/2","1","uud",1], 
     ["H"      ,"1837.47(-0.29)(0.20)","1837.47","-0.29","0.20"                ," "                    ,"0","","","H",-5.5],
     ["Neutron","1838.68366200(74)" ,"1838.68366200", "-0.00000074", "0.00000074","878.4(5)"           ,"0","1/2","1","udd",-5],  
     ["Eta`"   ,"1874.32(11)"       ,"1874.32","-0.11","0.11"                   ,"3.32(15) e -21"      ,"0","0","-","$u\overline{u}+\overline{d}d+s\overline{s}$",-22], 
     ["Phi"    ,"1995.035(31)"      ,"1995.035","-0.031","0.031"                ,"1.55(0,01) e -22"    ,"0","1","-","$s\overline{s}(most)$",-11],   
     ["c"      ,"2485(-39)(39)"     ,"2485","-39","39"                          ," "                       ,"+2/3","","","c",0],
     ["D"      ,"3670.4829677(11)"  ,"1837.47","-0.0000011","0.0000011"         ," "                    ,"0","","","D",-5.5],
     ["Tau"    ,"3477.23(23)"       ,"3477.23","-0.23","0.23"                   ,"290.3(5) e -15"           ,"-1","1/2","","tau",0],    
     ["b"      ,"8186(14)"          ,"8186","-14","14"                          ," "                       ,"-1/3","","","b",0],
     ["Higgs"  ,"244830(210)"       ,"244830","-210","210"                      ," "                       ,"0","0","","Higgs",0],
     ["t"      ,"337710(570)"       ,"337710","-570","570"                      ," "                       ,"+2/3","","","t",0]]

F=["#FFFFFF","#000000","#F60000","#05FB4F","#CFCF00","#000000","#07FCE4",
   "#F700D2","#00F73E","#7BB91F","#A9BF06","#047619","#047619","#789E20",
   "#CFCF00","#CF00B7","#EC61A9","#FA9805","#4200F6","#495999",
   "#B91F50","#CB4088","#F90404","#4D8E2F","#000000","#FFFFFF","#F50606",]

def Energie(i4,i3,i2,i1,i0,i_1,C):
    g[2][4]=i4; g[2][3]=i3; g[2][2]=i2; g[1][1]=i1; g[1][0]=i0; g[1][-1]=i_1
    E[0]=0; E[1]=0; E[2]=0; E[3]=0; E[4]=0; E[5]=0; E[6]=0; E[7]=0
    
    E_C_pos = -pi + 2*pi**(-1) - pi**(-3) + 2*pi**(-5) - pi**(-7) + pi**(-9) - pi**(-12) - 2*pi**(-14) 
    E_C_neg = 2*pi - pi**(-1) + E_C_pos    ;                             E[0]= pi**(-12) + 2*pi**(-14)                                      
    if C > 0:  E[0]= C * E_C_pos
    if C < 0:  E[0]= -C * E_C_neg
    
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

for i5 in [0]:    # for speed v  and Celestial bodies             
    for i4 in range(-2*O0,2*O0+1):                   # Select range for more particle: in range(-4,5):   
        for i3 in range(-2*O1,2*O1+1):               #            for Phi, Eta, tau, c    
           for i2 in range(-2*O1,2*O1+1):            # für Phi, Eta, ausreichend
                print("i4",i4,"i3",i3,"i2",i2,"i1",i1,"i_T1", i_T1)
                for i1 in range(-6,7):               # range(-6,7)  ist erforderlich      
                    for i0 in range(-6,7):  
                        for i_1 in range(-6,7):   
                            for C in range (-2,3):     
                                Energie(i4/2, i3/2, i2/2, i1/2, i0/2, i_1/2 ,C/2) 
                                if E[0] < 0: continue
                                if not Op == 5 and E[0] > 2000: continue     # Select any energy range.                  
                                if Op == 4 and (E[0]<1836 or E[0]>1839) : continue  # For H, Proton und Neutron: or E[0]<1836 or E[0]>1839
                                i_T += 1; flag = 0; N_T[j][1] = i_T
                                if N_T[j][0] == 0 : N_T[j][0] = i_T   
                                for j in range(1,26): 
                                    min_ = float(Obj[j][3]);  max_ = float(Obj[j][4])
                                    if (E[0]-float(Obj[j][2])<=1.0*max_) and (E[0]-float(Obj[j][2])>=1.0*min_):
                                        m =int(256+32*i4+4*i3+i2); i_T1 += 1;
                                        if Emax[j,m] <= E[0] :
                                            Emax[j,m]= E[0]; i_Emax[j,m]= i_T; Dmax[j,m,0]=i4; Dmax[j,m,1]=i3; Dmax[j,m,2]=i2
                                            Dmax[j,m,3]=i1; Dmax[j,m,4]=i0; Dmax[j,m,5]=i_1; Dmax[j,m,6]=C   
                                        if Emin[j,m] >= E[0] or Emin[j, m]== 0:
                                            Emin[j,m]= E[0]; i_Emin[j,m] = i_T; Dmin[j,m,0]=i4; Dmin[j,m,1]=i3; Dmin[j,m,2]=i2
                                            Dmin[j,m,3]=i1; Dmin[j,m,4]=i0; Dmin[j,m,5]=i_1; Dmin[j,m,6]=C 
                                        plt.plot(i_T,E[0], color=F[j], marker='.', markerfacecolor=F[j]); flag = 1
                                if E[0] > 0 and flag ==0:
                                    plt.plot([i_T,i_T+1], [E[0], E[0]], color="#C0BCBC")
print("..............   wait for Plot, several minutes   .................")
print("")
print("possible ET: ", i_T , "real ET: ",  i_T1)

print('{0:10}{1:5}{2:20}{3:16}{4:8}{5:5}{6:5}{7:5}{8:5}{9:5}{10:5}{11:5}'.format("particle","","mass Half-life",\
                    "     theory: E","   total ","  i4","  i3","  i2","  i1","  i0"," i-1","  C"), file=f)

for j in range (1,22):         #    bis   E<2000         [1,2,3,4,5,6,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26]:
    m_min = float(Obj[j][2]) + float(Obj[j][3]); m_max = float(Obj[j][2]) + float(Obj[j][4])
    p_g = len(Obj[j][2]); m_min =float(str(m_min)[:p_g]); m_max =float(str(m_max)[:p_g]);                flag = 0

    for m in range(1,512): 
        if i_Emax[j, m]==0: continue       #   Achtung  for i_1 in [0]:   (13**3) in (13**2) ersetzen
        E_mean = (Emax[j,m]+Emin[j,m])/2; Di_E=(i_Emax[j,m]-i_Emin[j,m]+1); i_Emax[j,0] +=abs(Di_E); D_i_E[j]=float(Di_E*100)/(5*13**3); D_i_N[j]=float(i_Emax[j,0]*100)/N_T[j][1]; i_Emax[j,-3]=float(Di_E*100)/(5**4*13**3)
        Emax[j,m]= float(str(Emax[j,m])[:p_g]); Emin[j,m]= float(str(Emin[j,m])[:p_g]); E_mean= float(str(E_mean)[:p_g])
        print("", file=f);   flag += 1
        print('{0:10}{1:5}{2:20}{3:14}{4:10}{5:11}{6:10}{7:2}{8:11}{9:10}{10:2}'.format(Obj[j][0],"",Obj[j][5],"    Σ abs(i_m)",i_Emax[j,0],\
                                  " ∆i_E/O_2 ","%.2f" % D_i_E[j]," %"," ∆i_E/N_T ","%.4f" % D_i_N[j]," % "),file=f);       
        print('{0:10}{1:5}{2:20}{3:16}{4:8}{5:5}{6:5}{7:5}{8:5}{9:5}{10:5}{11:5}'.format("","max   ", m_max,\
            Emax[j,m],i_Emax[j,m],Dmax[j,m,0]/2,Dmax[j,m,1]/2,Dmax[j,m,2]/2,Dmax[j,m,3]/2,Dmax[j,m,4]/2,Dmax[j,m,5]/2,Dmax[j,m,6]/2), file=f);   
        print('{0:10}{1:5}{2:20}{3:16}{4:8}{5:9}{6:7}{7:2}{8:12}{9:7}{10:2}'.format("","mean  ", Obj[j][1],E_mean,Di_E,"","","","","",""), file=f)                                                    
        print('{0:10}{1:5}{2:20}{3:16}{4:8}{5:5}{6:5}{7:5}{8:5}{9:5}{10:5}{11:5}'.format("","min   ", m_min,\
            Emin[j,m],i_Emin[j,m],Dmin[j,m,0]/2,Dmin[j,m,1]/2,Dmin[j,m,2]/2,Dmin[j,m,3]/2,Dmin[j,m,4]/2,Dmin[j,m,5]/2,Dmin[j,m,6]/2), file=f)
        if flag == 1:         
            if Op==1: X = [0,2,4,5, 3.5,4, 4,9, -13,-7,2, -11, 2,-14,-7, -11,-7, -12,-8,-3, -6,-9,-5.5,-22,-11,0,0,0]; Y=-20 ;  fs=12  #  text in plot
            if Op==2: X = [0,5,10,15, 20,5, 9,30, -35,-20,7, -40, -40,-23,5, -40,-22,  -65,-55,-45,  -35,-30, -20,-20, -40, 0,0,0,0]; Y=-20 ;  fs=12
            if Op==3: X = [0,1.2,2.5,3, 1,1, 0.7,1 ] ; Y = -4; fs=16   #,
            if Op==4: X = [0,2,4,5, 3,4, 4,8, -13,-7,2, -11, 2,-14,-7, -11,-7,  -12,-8,-3,  -0,-0,-0,-22,-11,0,0,0]; Y = -0; fs=16   #,
            if Op==5: X = [0,5,10,15, 20,5, 9,30, -35,-20,7, -40, -40,-23,5, -40,-22,  -65,-55,-45,  -35,-30, -20,-20, -40, 0,0,0,0]; X[j] *=2; Y = -0; fs=16   #,
            plt.text(i_Emin[j,m] + 10000 * X[j], E_mean + Y, Obj[j][9],fontsize=fs, color=F[j]); flag = 2 
            print(i_Emin[j,m] , E_mean, i_Emax[j,m])
    if flag==0: print('{0:10}{1:5}{2:20}{3:16}{4:8}'.format(Obj[j][0],"mean  ", Obj[j][1], "  ", " only with i4 > 1"), file=f);

if Op==1 or Op== 2:       #  for u,d,s  i4 <= 1   E<2000
    x_a= 0; x_m= i_T*1/5;  plt.ylabel('Energy in $m_e$'); plt.xlabel('N');  plt.xlim(-10000, i_T+30000); 
    i2 = float(2*pi)**2; i3 = float(2*pi)**3; i4 = float(2*pi)**4; i5 = 1/2*(i4+i3+i2); i6 = i4+i3+i2;  
    plt.plot([x_a,x_m],[i4,i4],'k',linewidth=1);plt.text(x_a,i4+15,'$(2\pi)^4$', fontsize=12, color='blue')
    plt.plot([x_a,x_m],[i3,i3],'k',linewidth=1);plt.text(x_a,i3+15,'$(2\pi)^3$', fontsize=12, color='blue')
    plt.plot([x_a,x_m],[i2,i2],'k',linewidth=1);plt.text(x_a,i2+15,'$(2\pi)^2$', fontsize=12, color='blue')
    plt.plot([x_a,x_m],[i5,i5],'k',linewidth=1);plt.text(x_a,i5+15,'$1/2((2\pi)^4+(2\pi)^3+(2\pi)^2)$', fontsize=12,color='blue')
    plt.plot([x_a,x_m],[i6,i6],'k',linewidth=1);plt.text(x_a,i6+15,'$(2\pi)^4+(2\pi)^3+(2\pi)^2$', fontsize=12,color='blue')

if Op==1:                               # Legend for u,d,s  i4 <= 1   E<2000              
    x_a = i_T *0.70; dx = i_T*0.08; i = -50               
    for j in [1,2,3,4,5,6,7,8,9,10,15,16,17,18,19]:        
        particle = str(Obj[j][9]) 
        # if j == 1  or j == 4  or j == 28:  particle=""        
        plt.text(x_a, i, Obj[j][0]); plt.text(x_a+dx,i,particle,color=F[j]); plt.text(x_a+2*dx, i, i_Emax[j,0]);
        D = "%.5f" % D_i_N[j] +" %"; plt.text(x_a+3*dx, i, D)
        i+= 100
    plt.text(x_a+2*dx, i, "  i  "); plt.text(x_a+3*dx, i, " ∆i/i_N ")    # plt.text(x_a+4*dx, i,"∆i/i_O1 ");

if Op==2:                                         # Legend for  i4 = 2    
    x_a = i_T *0.65; dx = i_T*0.08; i = -50                 
    for j in [4,5,6,7,8,9,10,15,16,17,18,19,20,21]:      
        particle = str(Obj[j][9]) 
        # if j == 1  or j == 4  or j == 28:  particle=""        
        plt.text(x_a, i, Obj[j][0]); plt.text(x_a+dx,i,particle,color=F[j]); plt.text(x_a+2*dx, i, i_Emax[j,0]);
        D = "%.5f" % D_i_N[j] +" %"; plt.text(x_a+3*dx, i, D)
        i+= 100
    plt.text(x_a+2*dx, i, "  i  "); plt.text(x_a+3*dx, i, " ∆i/i_N ")      
print(Op) 

if Op==3:                               # Legend for u,d,s             
    x_a= 0; x_m= i_T*1/5;  plt.ylabel('Energy in $m_e$'); plt.xlabel('N');  plt.xlim(-10000, i_T+30000); 
    i2 = float(2*pi)**2; i3 = float(2*pi)**3; i4 = float(2*pi)**4; i5 = 1/2*(i4+i3+i2); i6 = i4+i3+i2;  
    plt.plot([x_a,x_m],[i3,i3],'k',linewidth=1);plt.text(x_a,i3+15,'$(2\pi)^3$', fontsize=14, color='blue')
    plt.plot([x_a,x_m],[i2,i2],'k',linewidth=1);plt.text(x_a,i2+15,'$(2\pi)^2$', fontsize=14, color='blue')
    x_a = i_T *0.70; dx = i_T*0.10; i = 0               
    for j in [1,2,3,4,5,6,7]:        
        particle = str(Obj[j][9]) 
        # if j == 1  or j == 4  or j == 28:  particle=""        
        plt.text(x_a, i, Obj[j][0]); plt.text(x_a+dx,i,particle,color=F[j]); plt.text(x_a+2*dx, i, i_Emax[j,0]);
        D = "%.5f" % D_i_N[j] +" %"; plt.text(x_a+3*dx, i, D)       # D = "%.2f" % D_i_E[j] +" %"; plt.text(x_a+4*dx, i, D)
        i+= 20
    plt.text(x_a+2*dx, i, "  i  "); plt.text(x_a+3*dx, i, " ∆i/i_N ")    # plt.text(x_a+4*dx, i,"∆i/i_O1 "); 
print(Op) 
if Op==4:                    # Proton H-Atom Proton            
    x_a= 0; x_m= i_T*1/5;  plt.ylabel('Energy in $m_e$'); plt.xlabel('N');  plt.xlim(0, i_T); fs = 14
    plt.plot([1,i_T],[1836.1526734+1,1836.1526734+1],'k',linewidth=1); plt.text(1,1837,'$m_{Proton} + m_e$', fontsize=fs ,color='blue'); 
    x_a = i_T *0.70; dx = i_T*0.10; i = 0               
    for j in [17,19]:        
        plt.text(1050, float(Obj[j][2]), Obj[j][0],fontsize=fs, color=F[j]);
        particle = str(Obj[j][9]) 
        plt.text(x_a, i, Obj[j][0]); plt.text(x_a+dx,i,particle,color=F[j]); plt.text(x_a+2*dx, i, i_Emax[j,0]);
        D = "%.5f" % D_i_N[j] +" %"; plt.text(x_a+3*dx, i, D) ; D = "%.2f" % D_i_E[j] +" %"; plt.text(x_a+4*dx, i, D)
        i+= 20
print(Op) 
if Op==5:    # c   Deuterium    tau  
    x_a = i_T *0.65; dx = i_T*0.08; i = -50               
    for j in [1,5, 8,9,10,11,12,13,15,16,  17,18,19, 20,21,22, 23,24]:        
        particle = str(Obj[j][9]) 
        # if j == 1  or j == 4  or j == 28:  particle=""        
        plt.text(x_a, i, Obj[j][0]); plt.text(x_a+dx,i,particle,color=F[j]); plt.text(x_a+2*dx, i, i_Emax[j,0]);
        D = "%.5f" % D_i_N[j] +" %"; plt.text(x_a+3*dx, i, D)
        i+= 100
    plt.text(x_a+2*dx, i, "  i  "); plt.text(x_a+3*dx, i, " ∆i/i_N ")    # plt.text(x_a+4*dx, i,"∆i/i_O1 ");
print(Op) 

fig = plt.gcf()
fig.set_size_inches(10, 6)                     
fig.savefig(f_png, dpi=100)   # after several minutes
f.close()
plt.show()
exit()
