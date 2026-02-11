"""Entrypoint for Helmut :)"""

from polynome2pi.main import main

# =========================
# CONFIGURATION START
# =========================

ARGS = [
    "--sector", "222 u d s c"
    # "--no-show",
]

ARGSC = [
    "--Charge", "6",               #   Versuch  für  Argumente
]

# Examples:
#    minimal = "001 u d"
#    light = "011 u d s"
#    broad = "112 u d s nucleon"
#    nucleon = "111 H-atom"
#    heavy = "222 u d s c"
#    E112P = "112 E above 1700"

# =======================
# USER CONFIGURATION END
# =======================

if __name__ == "__main__":
    print("Running with arguments:", ARGS)
    main(ARGS)


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


#r_{apoapsis} = 696342 km \sqrt{1/2(2\pi)^5-1/2(2\pi)^4+(2\pi)^3} = 46006512 km $ \\
#$ Messung: \ 46.002 \ 10^6 km \ \ rel. Abweichung = 0.0001 $\\ \\
# $r_{periapsis} = 696342 km \sqrt{(2\pi)^5-0(2\pi)^4+(2\pi)^3} = 69775692 km $ \\
#$ Messung: \ 69.81 \ 10^6 km \ \ rel. Abweichung = 0.0005 $ \\ \\


#$r_{apoapsis} = 696342 km \sqrt{2(2\pi)^5+3(2\pi)^4-(2\pi)^3} = 107905705 km $\\
#$ Messung: \ 107.4128 \ 10^6 km \ \ rel. Abweichung = 0.004 $\\ \\
# $r_{periapsis} = 696342 km \sqrt{2(2\pi)^5+3(2\pi)^4+(2\pi)^3} = 109014662 km $\\
#$ Messung: \ 108.9088 \ 10^6 km \ \ rel. Abweichung = 0.001 $ \\ \\

# $r_{Venus} / r_{Mercury} = 6123.80 / 2448,57 = 2.50094$ & \hfill(3.4)\\

#1/2(2*pi)^5-1/2(2*pi)^4+(2*pi)^3
#(2*pi)^5-0(2*pi)^4+(2*pi)^3
   
#2(2*pi)^5+3(2*pi)^4-(2*pi)^3
#2(2*pi)^5+3(2*pi)^4+(2*pi)^3