# GAG Prospector
A Tool for Mining Glycosaminoglycan (GAG) Sequence Databases for Mass Spectrometry Data Analysis</br>
</br>
This software takes the numbers of monosaccharides and adduct ions as input parameters, and generates a database of GAG species within certain m/z and charge ranges.

## How to Use 
### Input parameters
<b>min_charge</b>: minimum charge</br>
<b>max_charge</b>: maximum charge</br>
<b>min_mz</b>: minimum m/z</br>
<b>max_mz</b>: maximum m/z</br>
<b>dHexA</b>: number of unsaturated hexenuronic acid</br>
<b>HexA</b>: number of hexenuronic acid</br>
<b>HexN</b>: number of hexosamine</br>
<b>HexNAc</b>: number of N-acetylhexosamine</br>
<b>Mann</b>: number of 2,5-anhydromannose</br>
<b>Ac</b>: maximum number of acetate</br>
<b>SO3</b>: maximum number of sulfate</br>
<b>NH4</b>: maximum number of ammonium</br>
<b>HOHloss</b>: water loss (maximum value: 1)</br>
<b>Na</b>: maximum number of sodium ion</br>
<b>K</b>: maximum number of potassium ion</br>
<b>Ca</b>: maximum number of calcium ion</br>
<b>Li</b>: maximum number of lithium ion</br>

### How to generate a GAG database
    min_charge = 4
    max_charge = 8
    min_mz = 200
    max_mz = 2000
    dHexA = 1
    HexA = 6
    HexN = 7
    HexNAc = 0
    Mann = 0
    Ac = 2
    SO3 = 20
    NH4 = 1
    HOHloss = 1
    exp = GagProspector(min_charge, max_charge, min_mz, max_mz, dHexA, HexA, HexN, HexNAc, Mann, Ac, SO3, NH4, HOHloss)
    db, _ = exp.build_database()
    if db:
        if len(db[0]) == 2:
            print "%s\t%s" % ("Mass", "dHexA HexA HexN HexNAc Mann Ac SO3 NH4 HOHloss Na K Ca Li")
            for mass, gag in db:
                print "%.4f\t%d %d %d %d %d %d %d %d %d %d %d %d %d" % (mass, gag.dHexA, gag.HexA, gag.HexN, gag.HexNAc, gag.Mann, 
                        gag.Ac, gag.SO3, gag.NH4, gag.HOHloss, gag.Na, gag.K, gag.Ca, gag.Li)
        else:
            print "%s\t%s\t%s" % ("m/z", "Charge", "dHexA HexA HexN HexNAc Mann Ac SO3 NH4 HOHloss Na K Ca Li")
            for mz, charge, gag in db:
                print "%.4f\t%d\t%d %d %d %d %d %d %d %d %d %d %d %d %d" % (mz, charge, gag.dHexA, gag.HexA, gag.HexN, gag.HexNAc,
                        gag.Mann, gag.Ac, gag.SO3, gag.NH4, gag.HOHloss, gag.Na, gag.K, gag.Ca, gag.Li)
    
