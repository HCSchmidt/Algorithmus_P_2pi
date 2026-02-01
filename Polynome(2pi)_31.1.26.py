import matplotlib.pyplot as plt  
import cmath
import matplotlib.colors as mcolors
import numpy as np
                    
Op = 5;    #  Selection of 5 plots with different resolutions and energy ranges. Option 1

if Op==1: J4=1; J3=1; J2=2; f_="Polynom 112 30.1.26.txt"; f_png = "Polynom 112 30.1.28.png" # slow (20 minutes)
if Op==2: J4=0; J3=1; J2=1; f_="Polynom 011 30.1.26.txt"; f_png = "Polynom 011 30.1.26.png" # fast, for u, d, s, pion, muon E<300
if Op==3: J4=0; J3=0; J2=1; f_="Polynom 001  30.1.26.txt"; f_png = 'Polynom 001 30.1.26.png' # fast, for u , d
if Op==4: J4=1; J3=1; J2=1; f_="Polynom 111 Atom 30.1.26.txt"; f_png = 'Polynom 111 Atom 30.1.26.png' # fast, for Proton H-Atom Neutron
if Op==5: J4=2; J3=2; J2=2; f_="Polynom 222 Tau 30.1.26.txt"; f_png = "Polynom 222 Tau 30.1.26.png" #  E<1500! slow, more particles with c,tau
f=open(f_,'w', encoding="utf8"); f=open(f_,'a+', encoding="utf8")
pi = cmath.pi; i_T = 0; i_T1 = 0; E = 0; i1=0;i0=0;i_1=0; mmax=0; m=0; E = [0]*10; g = [[0]*10]*10; j=0; D_i_N=[0]*300; 
Obj=[[0]*13]*35; N_E= [0]*100000; N_T= [[0]*100]*100; E_t=[0]*100000; X=[0]*35; D_i_c_=[0]*35; Cnt=[0]*520  
Emax = np.zeros((35, 520), dtype=float); Emin = np.zeros((35, 520), dtype=float)
i_Emax = np.zeros((35, 520), dtype=int); i_Emin = np.zeros((35, 520), dtype=int)
Dmax = np.zeros((35, 520, 7), dtype=int); Dmin = np.zeros((35, 520, 7), dtype=int)

Obj=[["Name"   ,"m_e"               ,"E",   "-SD",   "+SD"                      ,"Halbwertszeit T in sec" ,"Charge","Spin","P.","name"],
     ["e"      ,"1.00000000000(31)" ,"1.00000000000","-0.005","0.000"           ," "                       ,"-1","1/2","","e"],  
     ["u"      ,"4.18(-0.51)(0.96)" ,"4.18","-0.51","0.96"                      ," "                       ,"+2/3","","","u"],
     ["d"      ,"9.14(-0.33)(0.94)" ,"9.14","-0.33","0.94"                      ," "                       ,"-1/3","","","d"],
     ["s"      ,"182.8(-6.6)(16.8)" ,"182.8","-6.6","16.8"                      ," "                       ,"-1/3","","","s"],
     ["Muon"   ,"206.7682827(46)"   ,"206.7682827","-0.0000046","0.0000046"     ,"2.1969811(22) e -6"     ,"0","1","","muon"],
     ["Pion 0" ,"264.1430(9)"       ,"264.1430","-0.0009","0.0009"              ,"8.52(18) e -17"     ,"0","0","-","$u\overline{d}-\overline{u}d$"],        
     ["Pion +-","273.13243(35)"     ,"273.13243","-0.00035","0.00035"           ,"2.6033(5) e -8"     ,"+-1","0","-","$u\overline{u},\overline{d}d$"], 
     ["K +-"   ,"966.102(21)"       ,"966.102","-0.021","0.021"                 ,"1.2380(20) e -8"     ,"+-1","0","-","$u\overline{s},s\overline{u}$"], 
     ["KL 0"   ,"973.800(26)"       ,"973.800","-0.026","0.026"                 ,"5.116(21) e -8"      ,"0","0","-","$d\overline{s},s\overline{d}$  "],  
     ["KS 0"   ,"973.800(26)"       ,"973.800","-0.026","0.026"                 ,"8.954(4) e -11"      ,"0","","","$d\overline{s},s\overline{d}$"],   
     ["Eta"    ,"1072.139(35)"      ,"1072.139","-0.035","0.035"                ,"5 e -19"             ,"0","0","-","$u\overline{u}+\overline{d}d-2s\overline{s}$"], 
     ["Rho +-" ,"1506(1)"           ,"1506","-1","1"                            ,"4 e -24"             ,"-+1","1","-","$u\overline{u},\overline{d}d$"], 
     ["Rho 0"  ,"1517.14(49)"       ,"1517.14","-0.49","0.49"                   ,"4 e -24"             ,"0","1","-","$u\overline{u}-\overline{d}d$"], 
     ["Omega"  ,"1531.62(25)"       ,"1531.62","-0.25","0.25"                   ,"7.75(7) e -23"       ,"0","1","-","$u\overline{u}+\overline{d}d$"], 
     ["K* +-"  ,"1745.2(1)"         ,"1745.2","-0.1","0.1"                      ,"1.3 e -23"           ,"+-1","","","$d\overline{s},s\overline{d}$"],
     ["K* 0"   ,"1752.6(1)"         ,"1752.6","-0.1","0.1"                      ,"1.3 e -23"           ,"0","","","$d\overline{s},s\overline{d}$"],  
     ["Proton" ,"1836.152673426(32)","1836.152673426","-0.000000032","0.000000032"," "                 ,"1","1/2","1","uud"], 
     ["H"      ,"1837.47(-0.29)(0.20)","1837.47","-0.29","0.20"                ," "                    ,"0","","","H"],
     ["Neutron","1838.68366200(74)" ,"1838.68366200", "-0.00000074", "0.00000074","878.4(5)"           ,"0","1/2","1","udd"],  
     ["Eta`"   ,"1874.32(11)"       ,"1874.32","-0.11","0.11"                   ,"3.32(15) e -21"      ,"0","0","-","$u\overline{u}+\overline{d}d+s\overline{s}$"], 
     ["Phi"    ,"1995.035(31)"      ,"1995.035","-0.031","0.031"                ,"1.55(0,01) e -22"    ,"0","1","-","$s\overline{s}(most)$"],   
     ["c"      ,"2485(-39)(39)"     ,"2485","-39","39"                          ," "                       ,"+2/3","","","c",0],
     ["Tau"    ,"3477.23(23)"       ,"3477.23","-0.23","0.23"                   ,"290.3(5) e -15"           ,"-1","1/2","","tau"],    
     ["D 0"    ,"3649.38(10)"        ,"3649.38","-0.10","0.10"                  ,"4.101(15) e -13"     ,"+-1","0","-","$c\overline{u},u\overline{c}$"], 
     ["D +"    ,"3658.81(10)"        ,"3658.81","-0.10","0.10"                  ,"1.040(7) e -12"     ,"+-1","0","-","$c\overline{d},d\overline{c}$"], 
     ["Deuteron"  ,"3670.4829677(11)"  ,"3670.4829677","-0.0000011","0.0000011"    ," "                    ,"0","","","D",-5.5],
 
     ["DS +"  ,"3851.94(13)"        ,"3851.94","-0.13","0.13"                  ,"5.04(4) e -13"     ,"+-1","0","-","$c\overline{s},s\overline{c}$"], 
     ["Higgs"  ,"244830(210)"       ,"244830","-210","210"                      ," "                       ,"0","0","","Higgs"],
     ["t"      ,"337710(570)"       ,"337710","-570","570"                      ," "                       ,"+2/3","","","t"]]
     #  ab   223    oder   1222  erforderlich
     #["Eta_c"  ,"5839.3(9)"          ,"5839.3","-0.9","0.9"                     ,"2.06 e -23 "     ,"+-1","0","-","$c\overline{c},\overline{c}$"]],  
     #["b"      ,"8186(14)"           ,"8186","-14","14"                         ," "                       ,"-1/3","","","b"],
     #["Eta_b"  ,"18392.7(3.9)"       ,"18392.7","-3.9","3.9"                    ,""     ,"+-1","0","-","$b\overline{b},\overline{b}b"]

F=["#FFFFFF","#000000","#F60000","#05FB4F","#CFCF00","#000000","#07FCE4",
   "#F700D2","#00F73E","#7BB91F","#A9BF06","#047619","#047619","#789E20",
   "#CFCF00","#CF00B7","#EC61A9","#FA9805","#4200F6","#495999",
   "#B91F50","#CB4088","#F90404","#000000","#4D8E2F","#499999","#F50606","#146108"]

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

for i5 in [0]:      # for speed v  c Mesonen and Celestial bodies             
    for i4 in range(-2*J4,2*J4+1):                   # Select range for more particle: in range(-4,5):   
        for i3 in range(-2*J3,2*J3+1):               #            for Phi, Eta, tau, c    
           for i2 in range(-2*J2,2*J2+1):            
                print("i4",i4,"i3",i3,"i2",i2,"i1",i1,"i_T1", i_T1)
                for i1 in range(-6,7):               # range(-6,7) is required     
                    for i0 in range(-6,7):  
                        for i_1 in range(-6,7):   
                            for C in  range(-2,3): 
                                Energie(i4/2, i3/2, i2/2, i1/2, i0/2, i_1/2 ,C/2) 
                                if E[0] < 0: continue
                                m =int(256+32*i4+4*i3+i2);         # Counts in 2 pi with positive energy
                                if m > mmax : mmax =m; ct = 0                           
                                ct +=1; Cnt[mmax]=ct 
                                if Op == 5 and E[0] < 1500: continue   #2500: continue  # 
                                if Op == 4 and (E[0] < 1836 or E[0] > 1839): continue  # For H, Proton und Neutron: or E[0]<1836 or E[0]>1839
                                i_T += 1; flag = 0;
                                for j in range(1,27):     # für <6 1;26
                                    min_ = float(Obj[j][3]);  max_ = float(Obj[j][4])
                                    if (E[0]-float(Obj[j][2])<=1.0*max_) and (E[0]-float(Obj[j][2])>=1.0*min_):
                                        i_T1 += 1
                                        if Emax[j,m] <= E[0] :
                                            Emax[j,m]= E[0]; i_Emax[j,m]= i_T; Dmax[j,m,0]=i4; Dmax[j,m,1]=i3; Dmax[j,m,2]=i2
                                            Dmax[j,m,3]=i1; Dmax[j,m,4]=i0; Dmax[j,m,5]=i_1; Dmax[j,m,6]=C   
                                        if Emin[j,m] >= E[0] or Emin[j, m]== 0:
                                            Emin[j,m]= E[0]; i_Emin[j,m] = i_T; Dmin[j,m,0]=i4; Dmin[j,m,1]=i3; Dmin[j,m,2]=i2
                                            Dmin[j,m,3]=i1; Dmin[j,m,4]=i0; Dmin[j,m,5]=i_1; Dmin[j,m,6]=C 
                                        plt.plot(i_T, E[0], color=F[j],marker='.', markerfacecolor=F[j]); flag = 1   # 
                                if E[0] > 0 and flag ==0:
                                    plt.plot([i_T,i_T+1], [E[0], E[0]], color="#C0BCBC")
print("..............   wait for Plot, several minutes   .................")
print("")
print("possible ET: ", i_T , "real ET: ",  i_T1)
print("......with options 1 and 5, the process can take up to 60 minutes!.....")
print('{0:10}{1:5}{2:20}{3:16}{4:8}{5:5}{6:5}{7:5}{8:5}{9:5}{10:5}{11:5}'.format("particle","","",\
                    "     theory: E","   total ","  i4","  i3","  i2","  i1","  i0"," i-1","  C"), file=f)
for j in range (1,29):      
    m_min = float(Obj[j][2]) + float(Obj[j][3]); m_max = float(Obj[j][2]) + float(Obj[j][4])
    p_g = len(Obj[j][2]); m_min =float(str(m_min)[:p_g]); m_max =float(str(m_max)[:p_g]);                flag = 0

    for m in range(1,512): 
        if i_Emax[j, m]==0: continue
        if flag == 0: print(Obj[j][0], file=f);   flag += 1        

        E_mean = (Emax[j,m]+Emin[j,m])/2; Di_E=(i_Emax[j,m]-i_Emin[j,m]+1); i_Emax[j,0] +=abs(Di_E); Cnt_= Cnt[m]
        D_i_c=round((abs(Di_E))*100/Cnt_,5); i_Emax[j,516] += Cnt_; D_i_N[j] = float(i_Emax[j,0])*100/Cnt[m]; 
        Emax[j,m]= float(str(Emax[j,m])[:p_g]); Emin[j,m]= float(str(Emin[j,m])[:p_g]); E_mean = float(str(E_mean)[:p_g])

        print('{0:10}{1:5}{2:20}{3:16}{4:8}{5:5}{6:5}{7:5}{8:5}{9:5}{10:5}{11:5}'.format("","max   ", m_max,\
            Emax[j,m],i_Emax[j,m],Dmax[j,m,0]/2,Dmax[j,m,1]/2,Dmax[j,m,2]/2,Dmax[j,m,3]/2,Dmax[j,m,4]/2,Dmax[j,m,5]/2,Dmax[j,m,6]/2), file=f);   
        print('{0:10}{1:5}{2:20}{3:16}{4:8}{5:9}{6:7}{7:2}{8:12}{9:7}{10:2}'.format("","mean  ", Obj[j][1],E_mean,Di_E,"","","","","",""), file=f)                                                    
        print('{0:10}{1:5}{2:20}{3:16}{4:8}{5:5}{6:5}{7:5}{8:5}{9:5}{10:5}{11:5}'.format("","min   ", m_min,\
            Emin[j,m],i_Emin[j,m],Dmin[j,m,0]/2,Dmin[j,m,1]/2,Dmin[j,m,2]/2,Dmin[j,m,3]/2,Dmin[j,m,4]/2,Dmin[j,m,5]/2,Dmin[j,m,6]/2), file=f)
        print('{0:10}{1:5}{2:20}{3:15}{4:8}{5:5}{6:8}{7:9}{8:10}{9:2}'.format("","","","         ∆ abs(i)",abs(Di_E),\
                         "  Cts",Cnt_,"  ∆i/(2pi)", D_i_c," %"),file=f);        #"  %.4f" % D_i_N[0]

        if flag == 1:         
            if Op==1: X = [0,4,7,9, 6,9, 5,13, -21,-12,5, -19, 6,-28,-17, -21,-13, -20,-12,-10, -37,-9,-5.5,-22,-11,0,0,0]; Y=-20 ;  fs=12  #  text in plot
            if Op==2: X = [0,1.2,2.5,3, 1,1, 0.7,1]; Y = -4; fs=16  
            if Op==3: X = [0,0.2,0.2,0.2 ]; Y = -1; fs=16  
            if Op==4: X = [0,2,4,5, 3,4, 4,8, -13,-7,2, -11, 2,-14,-7, -11,-7,  -12,-8,-3,  -0,-0,-0,-22,-11,0,0,0]; Y = -0; fs=16   
            if Op==5: X = [0,5,10,15, 20,5, 9,30, -35,-20,7, -40, 25,15,4, 35,20,  -13,-10,-7,  20,20, 10,-10, -10, -20,0,0,0]; X[j] *=2; Y = -50; fs=12   #,
            plt.text(i_Emin[j,m] + 10000 * X[j], E_mean + Y, Obj[j][9],fontsize=fs, color=F[j]); flag = 2 
 
    if j > 1 and m == 511 and i_Emax[j,516] > 0:
        D_i_c = round( float(i_Emax[j,0])*100/i_Emax[j,516],5);   D = str(D_i_c);D_i_c_[j]="".rjust(10-len(D)," ") + D 
        print('{0:10}{1:5}{2:20}{3:15}{4:8}{5:5}{6:8}{7:9}{8:10}{9:2}'.format("","total","","             Σ ∆i",i_Emax[j,0],\
                           "  Cts",i_Emax[j,516], "  ∆i/(2pi)", D_i_c," %",'\n'),file=f);      

    if flag==0: print('{0:10}{1:5}{2:20}{3:16}{4:8}'.format(Obj[j][0],"mean  ", Obj[j][1], "  ", " only with i4 > 1"), file=f)

    x_a= 0; x_m= i_T*1/5;  plt.ylabel('Energy in $m_e$'); plt.xlabel('N');  plt.xlim(-10000, i_T+30000); 
    i2 = float(2*pi)**2; i3 = float(2*pi)**3; i4 = float(2*pi)**4; i5 = 1/2*(i4+i3+i2); i6 = i4+i3+i2; i7=5/2*i4-3/2*i3-1/2*i3; i8=2*i4+3*i3+2*i3; i9=2*i4+2*i3+2*i3; i10=3/2*i4+1/2*i3-1/2*i3;  
    if Op in [1] : plt.plot([x_a,x_m],[i4,i4],'k',linewidth=1); plt.text(x_a,i4+15,'$(2\pi)^4$', fontsize=12, color='blue')
    if Op in [5] : plt.plot([x_a,x_m],[i7,i7],'k',linewidth=1); plt.text(x_a,i7+15, '$5/2(2\pi)^4-3/2(2\pi)^3-1/2(2\pi)^2$', fontsize=12, color='blue')
    if Op in [5] : plt.plot([x_a,2*x_m],[i8,i8],'k',linewidth=1); plt.text(x_a,i8+15,'$2/2(2\pi)^4+3(2\pi)^3+3(2\pi)^2', fontsize=12, color='blue')
    if Op in [5] : plt.plot([x_a,3*x_m],[i9,i9],'k',linewidth=1); plt.text(x_a,i9+15,'$2(2\pi)^4+2(2\pi)^3+2(2\pi)^2$', fontsize=12, color='blue')
    if Op in [5] : plt.plot([x_a,3*x_m],[i10,i10],'k',linewidth=1); plt.text(x_a,i10+15,'$3/2(2\pi)^4+1/2(2\pi)^3-1/2(2\pi)^2$', fontsize=12, color='blue')
    if Op in [1,2] :plt.plot([x_a,x_m],[i3,i3],'k',linewidth=1); plt.text(x_a,i3+15,'$(2\pi)^3$', fontsize=12, color='blue')
    if Op in [1,2,3]: plt.plot([x_a,x_m],[i2,i2],'k',linewidth=1); plt.text(x_a,i2+15,'$(2\pi)^2$', fontsize=12, color='blue')
    if Op in [1]: plt.plot([x_a,x_m],[i5,i5],'k',linewidth=1); plt.text(x_a,i5+15,'$1/2((2\pi)^4+(2\pi)^3+(2\pi)^2)$', fontsize=12,color='blue')
    if Op in [1]: plt.plot([x_a,x_m],[i6,i6],'k',linewidth=1); plt.text(x_a,i6+15,'$(2\pi)^4+(2\pi)^3+(2\pi)^2$', fontsize=12,color='blue')

if Op==1:                               # Legend for u,d,s  i4 <= 1   E<2000              
    x_a = i_T *0.67; dx = i_T*0.08; i = -50               
    for j in [1,2,3,4,5,6,7,8,9,10,15,16,17,18,19]:        
        particle = str(Obj[j][9])
        plt.text(x_a, i, Obj[j][0]); plt.text(x_a+dx,i,particle,color=F[j]);
        plt.text(x_a+2*dx, i, i_Emax[j,0]); plt.text(x_a+3*dx, i, D_i_c_[j]) 
        i+= 100
    plt.text(x_a+2*dx, i, "  ∆i  "); plt.text(x_a+3*dx, i, " ∆i/(2pi) ")   

if Op==2:                               # Legend for u,d,s             
    x_m= i_T*1/5; plt.xlim(-10000, i_T+30000); x_a = i_T *0.70; dx = i_T*0.10; i = 0               
    for j in [1,2,3,4,5,6,7]:        
        particle = str(Obj[j][9]) 
        plt.text(x_a, i, Obj[j][0]); plt.text(x_a+dx,i,particle,color=F[j]);
        plt.text(x_a+2*dx, i, i_Emax[j,0]); plt.text(x_a+3*dx, i, D_i_c_[j])      
        i+= 20
    plt.text(x_a+2*dx, i, "  ∆i  "); plt.text(x_a+3*dx, i, " ∆i/(2pi) ")   

if Op==3:                               # Legend for u,d,s             
    x_m= i_T*1/5; plt.xlim(-1000, i_T+10000); x_a = i_T *0.89; dx = i_T*0.10; i = 0               
    for j in [1,2,3,4,5,6,7]:        
        particle = str(Obj[j][9]) 
        plt.text(x_a, i, Obj[j][0]); plt.text(x_a+dx,i,particle,color=F[j])
        plt.text(x_a+2*dx, i, i_Emax[j,0]); plt.text(x_a+3*dx, i, D_i_c_[j])      
        i+= 5
    plt.text(x_a+2*dx, i, "  ∆i  "); plt.text(x_a+3*dx, i, " ∆i/(2pi) ")   


if Op==4:                    # Proton H-Atom Proton            
    x_m= i_T*1/5;  plt.xlim(0, i_T); fs = 14; x_a = i_T *0.70; dx = i_T*0.10; i = 0               
    m_H = 1836.152673426+1
    plt.plot([1,i_T],[m_H,m_H],'k',linewidth=1); plt.text(1,1837,'$m_{Proton} + m_e$', fontsize=fs ,color='blue'); 
    for j in [17,19]:        
        plt.text(1050, float(Obj[j][2]), Obj[j][0],fontsize=fs, color=F[j]);
        particle = str(Obj[j][9]) 
        plt.text(x_a, i, Obj[j][0]); plt.text(x_a+dx,i,particle,color=F[j]);
        plt.text(x_a+2*dx, i, i_Emax[j,0]); plt.text(x_a+3*dx, i, D_i_c_[j])   
        i+= 20

if Op ==5 :    # c     tau  
    x_a = i_T *0.65; dx = i_T*0.085; i = 1500
    for j in [15,16,  17,18,19, 20,21,22, 23,24,25]:  #,26      
        particle = str(Obj[j][9]) 
        if j == 26:  particle=""        
        plt.text(x_a, i, Obj[j][0]); plt.text(x_a+dx,i,particle,color=F[j]);
        plt.text(x_a+2*dx, i, i_Emax[j,0]); plt.text(x_a+3*dx, i, D_i_c_[j])   
        i+= 100
    plt.text(x_a+2*dx, i, "  ∆i  "); plt.text(x_a+3*dx, i, " ∆i/(2pi) ")  

f.close()
fig = plt.gcf()
fig.set_size_inches(10, 6)                     
fig.savefig(f_png, dpi=100)   # after several minutes
plt.show()
exit()

#  https://github.com/users/HCSchmidt/projects/1
#  https://github.com/HCSchmidt/Algorithmus_P_2pi

 #  plt.text(i_Emin[j,m] , E_mean-20," ".ljust(X[j]) + Obj[j][9]+ " ".rjust(-X[j]),color=F[j],fontsize=12); flag = 2 

#j=10; X = [0,4,8,5, 3,4, 3,8, -13,-7,2, -10, 2,-14,-7,  -10,-6,  -12,-8,-3,  -6,-9,-5.5,-22,-11,0,0,0,0,0,0]  
#X[j] = 28
#print("WE"+"sdfg".center(int(15))+"fg")
#print(" ".ljust(X[j]) + Obj[j][9]+ " ".rjust(X[j]))
#
# for j in range(1,27): X[j]= int(float(Obj[j][2])); 
# X_E =  sorted(X)
# exit()
# import pandas as pd
# dat=0; df = pd.DataFrame({'particle':[],'':[],'meassure: m':[],'theory: E':[],'total':[],'i4':[],"i3":[],"i2":[],"i1":[],"i0":[],"i-1":[],"C":[]})
#        dat +=1; df.loc[dat] = [Obj[j][0]," max",m_max,Emax[j,m],i_Emax[j,m],Dmax[j,m,0]/2,Dmax[j,m,1]/2,Dmax[j,m,2]/2,Dmax[j,m,3]/2,Dmax[j,m,4]/2,Dmax[j,m,5]/2,Dmax[j,m,6]/2]
#        dat+=1;df.loc[dat] = ["","mean", Obj[j][1],E_mean,i_DeltaE,"i_DeltaE", i_DE , "‰", "max i_DeltaE", round(i_DE/(6**2),3),"‰",""]
#        dat +=1; df.loc[dat] = [""," min", m_min, Emin[j,m],i_Emin[j,m],Dmin[j,m,0],Dmin[j,m,1],Dmin[j,m,2],Dmin[j,m,3],Dmin[j,m,4],Dmin[j,m,5],Dmin[j,m,6]]

# df.round({"i4": 1, "i3": 1, "i2": 1, "i1": 1,"i0": 1, "i_1": 1, "C": 1})
# df['total'] = df['total'].apply(np.int32)
# # print(df)
# latex_code = df.to_latex(caption="Meine Tabelle aus Python",   label="tab:python_table") #,, column_format="lcccccccccccc"
# with open('mytable.tex','w') as tf:  tf.write(latex_code)
Obj=[["Name"   ,"m_e"               ,"","",""                               ,"MeV"              ,"Radius in fm","Charge","Spin","P.","Iso.","Iz","Quark","Marker"],
     ["e"      ,"1.00(1)"           ,"1.00","-0.01","0.01"                  ,"0.51099895069(16)",""                       ,""            ,"-1","1/2","","","","",'k.'], 
     ["u"      ,"4.18(-0.51)(0.96)" ,"4.18","-0.51","0.96"                  ,"2.16(-0.26)(0.49)",""                       ,""            ,"+2/3","","","","+1/2","u","r."],  
	 ["d"      ,"9.14(-0.33)(0.94)" ,"9.14","-0.33","0.94"                  ,"4.67(-0.17)(0.48)",""                       ,""            ,"-1/3","","","","-1/2","d","g."],  
	 ["s"      ,"182.8(-6.6)(16.8)" ,"182.8","-6.6","16.8"                  ,"93.4(-3.4)(8.6)"  ,""                       ,""            ,"-1/3","","","","-1","s","c."],    
	 ["c"      ,"2485(-39)(39)"     ,"2485","-39","39"                      ,"1270(-20)(20)"    ,""                       ,""            ,"+2/3","","","","1","c",""],
     ["b"      ,"8186(14)"          ,"8186","-14","14"                      ,"4183(7)"          ,""                       ,""            ,"-1/3","","","","-1","b",""],
     ["t"      ,"337710(570)"       ,"337710","-570","570"                  ,"172570(290)"      ,""                       ,""            ,"+2/3","","","","1","t",""],
     ["Higgs"  ,"244830(210)"       ,"244830","-210","210"                  ,"125110(110)"      ,""                       ,""            ,"0","0","","","","Boson",""],
     ["Muon"   ,"206.7682827(46)"   ,"206.7682827","-0.0000046","0.0000046" ,""                 ,"2.1969811(22) e -6" ,"R"           ,"-1","1/2","","","","","k."],  
     ["Pion 0" ,"264.1430(9)"       ,"264.1430","-0.0009","0.0009"          ,"134.9768(5)"      ,"8.52(18) e -17"     ,"R"           ,"0","0","-","1","0","uu-dd","y."],  
     ["Pion +-","273.13243(35)"     ,"273.13243","-0.00035","0.00035"       ,"139.57039(18)"    ,"2.6033(5) e -8"     ,"0,659(4)"    ,"+-1","0","-","1","+1","ud, ud","y."], 

     ["Eta "   ,"547,862(17)"       ,"547,862(17)","-0.021","0.021"          ,"547,862(17)"       ,"5.0(3) e -19"     ,"+-1","0","-","$u\overline{u}+\overline{d}d-2s\overline{s} /6$",-8], 
     ["Eta` "  ,"957,78(6)"         ,"957,78(6)","-0.021","0.021"             ,"957,78(6)"         ,"3.31(15) e -21"     ,"+-1","0","-","$u\overline{u}+\overline{d}d+2s\overline{s} /6$",-8], 
     ["Eta_c"  ,"5839.3(9)"         ,"5839.3","-0.9","0.9"             ,"2983,9(5)"         ,"2.06 e -23 "     ,"+-1","0","-","$c\overline{c},\overline{c}$",-8], 
     ["Eta_b"  ,"18392.7(3.9)"       ,"18392.7","-3.9","3.9"         ,"9398,7(2,0)"        ,""     ,"+-1","0","-","$b\overline{b},\overline{b}b",-8], 
     ["D +"   ,"3658.81(10)"        ,"3658.81","-0.10","0.10"           ,"1869,65(5)"      ,"1.040(7) e -12"     ,"+-1","0","-","$c\overline{d},d\overline{c}$",-8], 
     ["D 0"   ,"3649.38(10)"        ,"3649.38","-0.10","0.10"           ,"1864,83(5)"      ,"4.101(15) e -13"     ,"+-1","0","-","$c\overline{u},u\overline{c}$",-8], 
     ["D_S +" ,"3851.94(13)"        ,"3851.94","-0.13","0.13"           ,"1968,34(7)"      ,"5.04(4) e -13"     ,"+-1","0","-","$c\overline{s},s\overline{c}$",-8], 

     ["Eta"    ,"1072.139(35)"      ,"1072.139","-0.035","0.035"           ,"547.862(18)"      ,"5.0(3) e -19"            ,""            ,"0","0","-","0","0","uu+dd-2ss","m."], #u anti_s, s anti_u !!!! mittlere Lebensdauer  
     ["Eta`"    ,"1874.32(11)"       ,"1874.32","-0.11","0.11"               ,"957.78(6)"        ,"3.32(15) e -21"     ,""            ,"0","0","-","0","0","uu+dd+ss","m."], #u anti_s, s anti_u !!!! mittlere Lebensdauer   
     ["Rho +-" ,"1506(1)"           ,"506","-1","1"                        ,"770"              ,"4 e -24)"            ,""            ,"-+1","1","-","1","-+1","ud, ud","y."], #u anti_s, s anti_u !!!! mittlere Lebensdauer
     ["Rho 0"  ,"1517.14(49)"       ,"1517.14","-0.49","0.49"               ,"775.26(25)"       ,"4 e -24)"            ,""            ,"0","1","-","1","-","uu-dd","y."], #u anti_s, s anti_u !!!! mittlere Lebensdauer  
     ["Omega"  ,"1531.62(25)"       ,"1531.62","-0.25","0.25"               ,"782.66(13)"       ,"7.75(7) e -23"      ,""            ,"0","1","-","0","","uu+dd","y."], #u anti_s, s anti_u !!!! mittlere Lebensdauer  
     ["Phi"    ,"1995.035(31)"      ,"1995.035","-0.031","0.031"            ,"1019.461(0.016)"  ,"1.55(0,01) e -22"   ,""            ,"0","1","-","0","","ss(most)","m."], #u anti_s, s anti_u !!!! mittlere Lebensdauer  
     ["K +-"   ,"966.102(21)"       ,"966.102","-0.021","0.021"             ,"493.677(16)"      ,"1.2380(20) e -8"    ,"0.560(31)"   ,"+-1","0","-","1","+1/2","us, su","m."], #u anti_s, s anti_u !!!! mittlere Lebensdauer  
     ["KL 0"   ,"973.8900(26)"       ,"973.800","-0.026","0.026"             ,"497.611(13)"      ,"5.116(21) e -8"     ,"-0.077(10)"  ,"0","0","-","1/2","-1/2","ds, sd","m."], # fm^2; #d anti_s, s anti_d  s (KL)  
     ["KS 0"   ,"973.800(26)"       ,"973.800","-0.026","0.026"             ,"497.611(13)"      ,"8.954(4) e -11"     ,"-0.077(10)"  ,"0","","","","","ds, sd","m."], # fm^2; #d anti_s, s anti_d  s (KL)  
     ["K* +-"  ,"1745.2(1)"         ,"1745.2","-0.1","0.1"                  ,"891.8"            ,"1.3 e -23"          ,""            ,"+-1","","","","","us, su","m."],  #u anti_s, s anti_u#  
     ["K* 0"   ,"1752.6(1)"         ,"1752.6","-0.1","0.1"                  ,"895.6"            ,"1.3 e -23"          ,""            ,"0","","","","","ds, sd","m."],  #d anti_s, s anti_d#   
     ["Neutron","1838.68366200(74)" ,"1838.68366200", "-0.00000074", "0.00000074"    ,""        ,"878.4(5)"               ,""            ,"0","1/2","1","1/2","-1/2","udd","b."],  
     ["Proton" ,"1836.152673426(32)","1836.152673426","-0.000000032","0.000000032"   ,""        ,""                      ,"0.8409(4)"   ,"1","1/2","1","1/2","+1/2","uud","b."], 
     ["Tau"    ,"3477.23(23)"       ,"3477.23","-0.23","0.23"               ,""                 ,"290.3(5) e -15"     ,"R"           ,"-1","1/2","","","","","k."],    
     ["W +-"   ,"157278.6(26.0)"    ,"157278.6","-26.0","26.0"              ,"80369.2(133)"     ,""                       ,""            ,"+-1","","","","","",""],
     ["Z 0"    ,"178450.4(4.5)"     ,"178450.4","-4.5","4.5"                ,"91188.0(20)"      ,""                       ,""            ,"0","","","","","",""],
     ["H"      ,"1837.47(-0.29)(0.20)","1837.47","-0.29","0.20"             ,""                 ,""                       ,""            ,"0","","","","","","b."],
     ["He"     ,"7296.2971(36)"     ,"7296.2971","-0.0036","0.0036"         ,""                 ,""                       ,""            ,"0","","","","","",""],
     ["Earth"  ,""                  ,"","",""                               ,""                 ,"86400"                  ,"6378137.0 m" ,"6356752.314","C","","","","","",""],
     ["Moon"   ,""                  ,"","",""                               ,""                 ,"27.3*24*3600"           ,"3474.8/2 m"  ,"C","","","","","",""]]

print("i_Emax[j,m]", i_Emax[j,m] , "10000*float(Obj[j][10]" ,10000*float(Obj[j][10]),"E_mean-20", E_mean-20, "Obj[j][9]",Obj[j][9])
#import sys  #from datetime import datetime; DT=datetime.today().strftime('%H:%M:%S'); print(DT); print(DT, file=f)

import pandas as pd
dat=0
dat +=1; data[dat] = ["particle","","mass Half-life", "theory: E","   total ","  i4","  i3","  i2","  i1","  i0"," i-1","  C"]
dat +=1; data[dat] = ["","max   ", m_max, Emax[j,m],i_Emax[j,m],Dmax[j,m,0]/2,Dmax[j,m,1]/2,Dmax[j,m,2]/2,Dmax[j,m,3]/2,Dmax[j,m,4]/2,Dmax[j,m,5]/2,Dmax[j,m,6]/2],
dat +=1; data[dat] = ["","mean  ", Obj[j][1],E_mean,Di_E,"","","","","",""],
dat +=1; data[dat] = ["","min   ", m_min, Emin[j,m],i_Emin[j,m],Dmin[j,m,0]/2,Dmin[j,m,1]/2,Dmin[j,m,2]/2,Dmin[j,m,3]/2,Dmin[j,m,4]/2,Dmin[j,m,5]/2,Dmin[j,m,6]/2]
for i in range [1,dat]:
    df = pd.DataFrame(data)
df.to_excel('output.xlsx', index=False)
