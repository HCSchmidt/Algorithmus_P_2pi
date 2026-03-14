from polynome2pi.main import main
# =========================
# CONFIGURATION START
# =========================
#    Selection of the experiment
Exp = [ "000 u d", "012 u d s", "112 u d s nucleon",
      "111 H-atom", "112 E > 1700", "222 u d s c",
      "333 E > 1700", "333 E < 2000",  
    ]
ARGS = [
        "--sector", Exp[6]
    ]
# select Charge  [" ","1","2/3","1/3","0","-1/3","-2/3","-1","+-1"] 
ARGC = [
        "--Charge", " ",     
    ]
# select H Atom  ["H","N"] or [" "," "] to see the Proton
ARGH = [
        "--HAtom", "H", "N"   
    ]

if __name__ == "__main__":
        print("Running with arguments:", ARGS, ARGC, ARGH)
        main(ARGS, ARGC, ARGH)
