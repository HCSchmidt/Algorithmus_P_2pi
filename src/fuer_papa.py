"""Entrypoint for Helmut :)"""

from polynome2pi.main import main

# =========================
# CONFIGURATION START
# =========================

ARGS = [
    "--sector", "broad"
    # "--no-show",
]

ARGSC = [
    "--Charge", "6",               #   Versuch  für  Argumente
]

# Examples:
# ARGS = ["--sector", "minimal"]
# ARGS = ["--sector", "nucleon"]
# ARGS = ["--sector", "broad"]
# ARGS = ["--sector", "heavy", "--base-name", "tau_scan"]
# ARGS = ["--sector", "light"]

# =======================
# USER CONFIGURATION END
# =======================

if __name__ == "__main__":
    print("Running with arguments:", ARGS)
    main(ARGS)


#  Alex
#  ich würde vorschlagen die
#  die Sektoren umzubenennen amit es zum Paper passt 
#
#  minimal >> 001 u_d 
#  light   >> 011 u_d_s 
#  broad   >> 112 nucleon
#  heavy   >> 222 c_tau
#  nucleon >> 112 H-Atom 
#
#  In constants habe ich die Kategorie ergänzt, damit es übersichtlicher ist  
#
#  ich möchte die bitten auch eine text datei .txt dazuzufügen
#  wo die events in 3 zeilen mit mean, min und max übereinander liegen 



#  mir ist jetzt erst aufgefallen, dass die .csv  an manchen stellen z.B. pion 0 komischerweise
#  die Zahlen nicht als zahlen sonder als Text angezeigt werden.
#  in notebook++ ist alles korrekt. Es liegt wohl am openoffice mit dem ich .csv öffne. 
#  er macht aus 264.1430   264143   und offensichtlich weil er die hinten  die 0 nicht apzetier
#  
# es liegt an den 3 Zeilen um die Zahlen  mean -min und max auf die gleiche Länge zu bringen.
#
#  p_g = len(Obj[j][2])
#  m_min = float(str(m_min)[:p_g])
#  m_max = float(str(m_max)[:p_g])
#
# eingentlich ein einfaches Problem.
# Meine Frage deshab funktioniert dies Exel oder nicht??
#





# und super wäre eine neue Liste von Objekten 

AE = 0.387098   # astronomische Einheit
def build_sun_table():
    Obj_s = [
#  "Catégorie"
        ["Name"   ,"Radius in Km","polar Radius","Neigung"   ,"radius_Orbit"   ,"-SD_r_O"     ,"SD_r_O"       ,"period in day","Ekliptik in °" ],
#  "Sun" 
        ["Sun"    , 696342       ,    0        ,     0       , 0                 , 0            , 0             , 25.38         ,      0       ],
#  "Planet" 
        ["Mercury", 4881 /2      , 4876.6 /2   ,    0.034    ,  57900000        , 0.308 * AE   , 0.467 * AE    , 87.969        ,      7.004    ],
        ["Venus"  , 12103.6 /2   , 12103.6 /2  , 177.36      , 108200000        , 0,718 * AE   , 0,728 * AE    , 224.701       ,      3.395    ],
        ["Earth"  , 6378137.0    , 6356752.314 , 23.44       , 149600000        , 0.983 * AE   , 1.017 * AE    , 365.256       ,      0.0001   ],
        ["Mars"   , 6792.4 /2    , 6752.4 /2   , 25.19       , 227990000        , 1.382 * AE   , 1.666 * AE    , 686.980       ,   1.8506      ],
#  "Moon"     
        ["Moon"   , 3474 /2      , 3474 /2     , 6.68        , 384400           , 363300       , 405500        , 27.3217       ,   5.145       ],
    ]                    

# das ist nur der Einstieg um das Sonnensystem einzubinden. 
# Ob das zu machen ist, weiss ich nicht. Der Algorithmus müsste angepasst werden 
# wohl mit einer eigene engine mit i5 (2pi)^5 . vielleich in den nächsten Monaten
# ich habe das schon mal in visula basic gemacht mit einem Plot wo man die Planeten als video 
# laufen lassen kann. Die Standardabweichungen geben die minimale und maximale Entfernung an

