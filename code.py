import pandas as pd
while True:
    file_pathway = "D:\data\LambdaPy\demo_input.csv"
    try:
        data = pd.read_csv(file_pathway)
        print("File successfully loaded.")
        break
    except FileNotFoundError:
        print("\nError: File not found. Please check the pathway and try again.\n")
    except pd.errors.EmptyDataError:
        print("\nError: The file is empty.\n")
    except pd.errors.ParserError:
        print("\nError: Unable to parse the CSV file. Please ensure it is formatted correctly.\n")
    except Exception as e:
        print("\n An unexpected error occurred:", e, '\n')


while True:
    reaction_pathway = 1.1
    try:
        reaction_pathway = float(reaction_pathway)
        if reaction_pathway in [0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 3.1, 3.2, 3.3, 3.4, 4.1, 4.2, 5, 6]:
            break
        else:
            print("Please enter a valid value for an electron acceptor reaction pathway from the provided list")
    except ValueError:
        print("Invalid input. Please enter a number.")


# INPUT 3
while True:
    pH_Input = 7
    try:
        pH_Input = int(pH_Input)
        if 0 <= pH_Input <= 14:
            break
        else:
            print("Value must be between 0 and 14.")
    except ValueError:
        print("Invalid input. Please enter an integer value.")


import pandas as pd
import numpy as np
import math 

pd.set_option('display.max_columns', 200)

CHEMICAL_ELEMENTS = ["C", "H", "N", "O", "P", "S"]

pathway = file_pathway
ui = reaction_pathway

def get_compositions(df):
    chemical_compositions = None
    formulas = None
    if "C" in df.columns:
        tdf = df[df["C"] > 0]
        if "C13" in df.columns:
            tdf = tdf[tdf["C13"] == 0]
        chemical_compositions = np.array(tdf[CHEMICAL_ELEMENTS])
        formulas = tdf["MolForm"].tolist()

    else:
        raise ValueError("Either columns for compositions (e.g., C, H, N, ...) column is required.")
    if "Z" in df.columns:
        chemical_compositions = np.column_stack((chemical_compositions, df["Z"]))
    else:
        chemical_compositions = np.column_stack((chemical_compositions, np.zeros(len(df))))
    return {
        "chemical_compositions": chemical_compositions,
        "formulas": formulas
    }


Compounds = ['H2O',
 'HCO3-',
 'NH4+',
 'HPO4',
 'HS-',
 'H+',
 'e-',
 'O2',
 'NO3-',
 'NO2',
 'N2',
 'Fe3+',
 'Fe2+',
 'H2',
 'SO4',
 'H2S',
 'SO3',
 'S',
 'S2O3',
 'Fe(OH)3',
 'FeOOH',
 'Fe3O4',
 'CH4',
 'CH3COO-',
 'MnO2',
 'Mn2+',
 'Biomass']

GE_of_compounds = [-237.2,
 -586.9,
 -79.37,
 -1089.1,
 12.05,
 0.0,
 0.0,
 16.5,
 -111.3,
 -32.2,
 18.19,
 -4.6,
 -78.87,
 0.0,
 -744.63,
 -33.4,
 -486.6,
 0.0,
 522.5,
 -690.0,
 -489.8,
 -1012.6,
 -34.06,
 -392.0,
 -465.14,
 -228.0,
 -67.0]

SDG = pd.DataFrame({"Compounds" : Compounds, "Standard Gibbs Energy of this compound in KJ/mol" : GE_of_compounds})

# 0) Oxygen (Default set to O2 as e acceptor)
# **0 (Default). O2 + 4H+ + 4e- ----> 2H2O**


# 1) Nitrogen Compounds, Nitrates and Nitrites
# **1.1.  NO3- + 10H+ + 8e- ---> NH4+ + 3H2O %**

# **1.2.  NO3- + 2H+ + 2e- ---> NO2- + H2O**

# **1.3.  NO3- + 6H+ + 5e- ---> (1/2)N2 + 3H2O**

# **1.4   NO2- + 4H+ + 3e- ---> (1/2)N2 + 2H2O**

# **1.5.  NO2- + 8H+ + 6e- ---> NH4+ + 2H2O**

# **1.6.  N2 + 8H+ + 6e- ---> 2NH4+**


# 2) Sulphur compounds, Sulphates and Sulphites
# **2.1.  SO4^2- + (9/2)H+ + 8e- ---> (1/2)H2S + (1/2)HS- + 4H2O**

# **2.2 SO4^2- + 2H+ + 2e- ---> SO3^2- + H2O**

# **2.3. SO4^2- + 5H+ + 4e- ---> (1/2)S2O3^2- + (5/2)H2O**

# **2.4. SO4^2- + 8H+ + 6e- ---> S + 4H2O**

# **2.5. SO4^2- + 9H+ + 8e- --> HS- + 4H2O**

# **2.6. SO3^2- + (15/2)H+ + 6e- ---> (1/2)H2S + (1/2)HS- + 3H2O**


# 3) Iron compounds, ferrous and ferric
# **3.1. Fe(OH)3 + 3H+ + e- --> Fe2+ + 3H2O**

# **3.2. FeOOH + 3H+ + e- --> Fe2+ + 2H2O**

# **3.3. Fe3O4 + 8H+ + 2e- --> 3Fe2+ + 4H2O**

# **3.4.  Fe3+ + e- ---> Fe2+**


# 4) Bicarbonate and Hydrogen ion
# **4.1. HCO3- + 9H+ + 8e- --> CH4 + 3H2O**

# **4.2. H+ + e- ---> (1/2)H2**


# 5) Acetate
# **5  CH3COO- + 9H+ + 8e- --> 2CH4 + 2H2O**


# 6) Manganese
# **6 MnO2 + 4H+ + 2e- --> Mn2+ + 2H2O**



def getThermoStoich(chemForm):

    a = chemForm[0]
    b = chemForm[1]
    c = chemForm[2]
    d = chemForm[3]
    e = chemForm[4]
    f = chemForm[5]
    z = chemForm[6]

     # Step 1a) stoichD: stoichiometries for an electron donor
    ySource = -1
    yH2o = -(3*a+4*e-d)
    yHco3 = a
    yNh4 = c
    yHpo4 = e
    yHs = f
    yH = 5*a+b-4*c-2*d+7*e-f
    yE = -z+4*a+b-3*c-2*d+5*e-2*f
    stoichD = [ySource,yH2o,yHco3,yNh4,yHpo4,yHs,yH,yE]
    stoichD[8:28] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,0 , 0, 0, 0, 0,0]

     # Step 1b) stoichA: stoichiometries for an electron acceptor 
    stoichA = np.zeros(28)      
    if ui == 1.1:
        stoichA = np.zeros(28)
        stoichA[9] = -1
        stoichA[6] = -10
        stoichA[7] = -8
        stoichA[3] = 1
        stoichA[1] = 3

    elif ui == 1.5:
        stoichA = np.zeros(28)
        stoichA[10] = -1
        stoichA[6] = -8
        stoichA[7] = -6
        stoichA[3] = 1
        stoichA[1] = 2

    elif ui == 1.6:
        stoichA = np.zeros(28)
        stoichA[11] = -1
        stoichA[6] = -8
        stoichA[7] = -6
        stoichA[3] = 2
        stoichA[1] = 0

    elif ui == 3.4:
        stoichA = np.zeros(28)
        stoichA[12] = -1
        stoichA[6] = 0
        stoichA[7] = -1
        stoichA[13] = 1
        stoichA[1] = 0

    elif ui == 4.2:
        stoichA = np.zeros(28)
        stoichA[6] = -1
        stoichA[7] = -1
        stoichA[14] = 0.5


    elif ui == 1.2:
        stoichA = np.zeros(28)
        stoichA[9] = -1
        stoichA[6] = -2
        stoichA[7] = -2
        stoichA[10] = 1
        stoichA[1] = 1

    elif ui == 1.3:
        stoichA = np.zeros(28)
        stoichA[9] = -1
        stoichA[6] = -6
        stoichA[7] = -5
        stoichA[11] = 0.5
        stoichA[1] = 3

    elif ui == 1.4:
        stoichA = np.zeros(28)
        stoichA[10] = -1
        stoichA[6] = -4
        stoichA[7] = -3
        stoichA[11] = 0.5
        stoichA[1] = 2

    elif ui == 2.1:
        stoichA = np.zeros(28)
        stoichA[15] = -1
        stoichA[6] = -4.5
        stoichA[7] = -8
        stoichA[16] = 0.5
        stoichA[5] = 0.5
        stoichA[1] = 4

    elif ui == 2.6:
        stoichA = np.zeros(28)
        stoichA[17] = -1
        stoichA[6] = -7.5
        stoichA[7] = -6
        stoichA[16] = 0.5
        stoichA[5] = 0.5
        stoichA[1] = 3

    elif ui == 2.2:
        stoichA = np.zeros(28)
        stoichA[15] = -1
        stoichA[6] = -2
        stoichA[7] = -2
        stoichA[17] = 1
        stoichA[1] = 1

    elif ui == 2.4:
        stoichA = np.zeros(28)
        stoichA[15] = -1
        stoichA[6] = -8
        stoichA[7] = -6
        stoichA[18] = 1
        stoichA[1] = 4

    elif ui == 2.3:
        stoichA = np.zeros(28)
        stoichA[15] = -1
        stoichA[6] = -5
        stoichA[7] = -4
        stoichA[19] = 0.5
        stoichA[1] = 2.5


    elif ui == 0:
        stoichA = np.zeros(28)
        stoichA[8] = -1
        stoichA[6] = -4
        stoichA[7] = -4
        stoichA[1] = 2


    elif ui == 3.1:
        stoichA = np.zeros(28)
        stoichA[20] = -1
        stoichA[6] = -3
        stoichA[7] = -1
        stoichA[13] = 1
        stoichA[1] = 3

    elif ui == 3.2:
        stoichA = np.zeros(28)
        stoichA[21] = -1
        stoichA[6] = -3
        stoichA[7] = -1
        stoichA[13] = 1
        stoichA[1] = 2

    elif ui == 3.3:
        stoichA = np.zeros(28)
        stoichA[22] = -1
        stoichA[6] = -8
        stoichA[7] = -2
        stoichA[13] = 3
        stoichA[1] = 4

    elif ui == 2.5:
        stoichA = np.zeros(28)
        stoichA[15] = -1
        stoichA[6] = -9
        stoichA[7] = -8
        stoichA[5] = 1
        stoichA[1] = 4

    elif ui == 4.1:
        stoichA = np.zeros(28)
        stoichA[2] = -1
        stoichA[6] = -9
        stoichA[7] = -8
        stoichA[23] = 1
        stoichA[1] = 3

    elif ui == 5:
        stoichA = np.zeros(28)
        stoichA[24] = -1
        stoichA[6] = -9
        stoichA[7] = -8
        stoichA[23] = 2
        stoichA[1] = 2

    elif ui == 6:
        stoichA = np.zeros(28)
        stoichA[25] = -1
        stoichA[6] = -4
        stoichA[7] = -2
        stoichA[26] = 1
        stoichA[1] = 2

    # Step 1c) stoichCat: stoichiometries for catabolic reaciton 
    yEd = stoichD[7]
    yEa = stoichA[7]
    stoichCat = [(stoichD[i] - ((yEd/yEa)*stoichA[i])) for i in range(len(stoichD))]

    # Step 2a) stoichAnStar: stoichiometries for anabolic reaciton (N source = NH4+)
    chemFormBiom = [1, 1.8, 0.2, 0.5, 0, 0, 0]
    aB = chemFormBiom[0]
    bB = chemFormBiom[1]
    cB = chemFormBiom[2]
    dB = chemFormBiom[3]
    eB = chemFormBiom[4]
    fB = chemFormBiom[5]
    zB = chemFormBiom[6]
    ySource = -1
    yH2o = -(3*aB+4*eB-dB)
    yHco3 = aB
    yNh4 = cB
    yHpo4 = eB
    yHs = fB
    yH = 5*aB+bB-4*cB-2*dB+7*eB-fB
    yE = -zB+4*aB+bB-3*cB-2*dB+5*eB-2*fB
    stoichAnStarB = [ySource,yH2o,yHco3,yNh4,yHpo4,yHs,yH,yE]
    stoichAnStarB[8:28] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,0 , 0, 0, 0, 0,0] # add additional components: e-acceptor and biomass  
    stoichAnStarB = [-x for x in stoichAnStarB]
    stoichAnStarB[27] = stoichAnStarB[0]
    stoichAnStarB[0] = 0

     # Step 2b) "overall" anabolic reaction
    stoichAnStar = [(stoichAnStarB[i] + ((1/a)*stoichD[i])) for i in range(len(stoichAnStarB))]
    yEana = stoichAnStar[7]
    if yEana > 0:
        stoichAn = [(stoichAnStar[i] - ((yEana/yEa)*stoichA[i])) for i in range(len(stoichAnStar))]
    elif yEana < 0:
        stoichAn = [(stoichAnStar[i] - ((yEana/yEd)*stoichD[i])) for i in range(len(stoichAnStar))]
    else:
        stoichAn = stoichAnStar


     # Step 3: get lambda

    # - estimate delGd0 using LaRowe and Van Cappellen (2011)

    ne = -z+4*a+b-3*c-2*d+5*e-2*f  # number of electrons transferred in D 
    nosc = -ne/a+4  # nominal oxidataion state of carbon 
    delGcox0 = 60.3-28.5*nosc  # kJ/C-mol
    delGd0 = delGcox0*a*abs(stoichD[0])  # kJ/rxn


    # - estimate delGf0 for electron donor
    delGf0_D_zero = 0
    delGf0_zero = [delGf0_D_zero, -237.2, -586.9, -79.37, -1089.1, 12.05, 0, 0, 16.5, -111.3, -32.2, 18.19, -4.6, -78.87, 0, -744.63, -33.4, -486.6, 0, 522.5, -690, -489.8, -1012.6, -34.06, -392, -465.14, -228, -67] 
    delGcox0_zero = np.dot(delGf0_zero, stoichD)
    delGf0_D_est = (delGd0-delGcox0_zero)/stoichD[0]  
    delGf0 = delGf0_zero
    delGf0[0] = delGf0_D_est


    # - standard delG at pH=0
    delGcat0 = np.dot(delGf0, stoichCat)
    delGan0 = np.dot(delGf0, stoichAn)

    # - stadard delG at pH=7
    R = 0.008314
    T = 298.15
    iProton = 6
    pH_value = 10**(-pH_Input)
    delGd = delGd0 + (R*T*stoichD[iProton])*(math.log(pH_value))
    delGcox = delGd / a
    delGcat = delGcat0+ (R*T*stoichCat[iProton])*(math.log(pH_value))
    delGan = delGan0+ (R*T*stoichAn[iProton])*(math.log(pH_value))


    # The Thermodynamic Electron Equivalents Model (TEEM)
    eta = 0.43
    delGsyn = 200
    if math.isnan(delGan0) and math.isnan(delGan):
        lambda0 = math.nan
        lambda_ = math.nan
        stoichMet = [math.nan] * len(stoichCat)
        delGdis0 = math.nan
        delGdis = math.nan
    else:
        if delGan < 0:
            m = 1
        else:
            m = -1
        if delGan0 < 0:
            m0 = 1
        else:
            m0 = -1
        lambda0 = ((delGan0*(eta**m0))+delGsyn)/(-delGcat0*eta)
        lambda_ = (delGan*eta**m+delGsyn)/(-delGcat*eta)
        if lambda_ > 0:
            stoichMet = [lambda_*stoichCat[i] + stoichAn[i] for i in range(len(stoichCat))]
        else:
            stoichMet = stoichAn
        delGdis0 = np.dot(lambda0, -delGcat0) - delGan0
        delGdis = np.dot(lambda_, -delGcat) - delGan

    #Organizing all the calculated values in a dataframe
    CUE = stoichMet[-1] * 1 / (abs(stoichMet[0]) * a)
    NUE = stoichMet[-1] * 0.2 / (abs(stoichMet[0]) * c + abs(stoichMet[3]) * 1)
    TER = (NUE / CUE) * (1 / 0.2)

    names = [""] * 62
    names[0:12] = ["CUE", "NUE", "TER","delGcox0","delGd0","delGcat0","delGan0","delGdis0","lambda0",
               "delGcox","delGd","delGcat","delGan","delGdis","lambda"]
    stoich_colnames = ["donor","h2o","hco3","nh4","hpo4","hs","h","e","O2", "NO3","NO2", "N2", "Fe3+", "Fe2+", "H2", "SO4", "H2S", "SO3", "S", "S2O3", "Fe(OH)3", "FeOOH", "Fe3O4", "CH4", "CH3COO-", "MnO2", "Mn2+" , "biom"]
    stoich_types = ["stoichD","stoichA","stoichCat","stoichAn","stoichMet"]
    for i in range(len(stoich_types)):
        names[(i*28+15):(i*28+25)] = [stoich_types[i] + "_" + colname for colname in stoich_colnames]
    all_values = [CUE, NUE, TER, delGcox0, delGd0, delGcat0, delGan0, delGdis0, lambda0, delGcox, delGd, delGcat, delGan, delGdis, lambda_]
    all_values.extend(stoichD)
    all_values.extend(stoichA)
    all_values.extend(stoichCat)
    all_values.extend(stoichAn)
    all_values.extend(stoichMet)
    df1 = pd.DataFrame(data = [all_values], columns = names)
    return df1


def get_lambda(compositions):
    b1 = pd.DataFrame() 
    for i in compositions:
        b = getThermoStoich(i)  
        b = pd.concat([b1, b])
        b1 = b

    return b1





try:
    gt = get_compositions(pd.read_csv(pathway))
except:
    print("Please enter a correct pathway for a csv file")
    

try:
    df_lambda = get_lambda(gt["chemical_compositions"])
    df_lambda.to_csv(f"Output_with_Rxn_no_{ui}.csv", index=False)
except:
    print(f"Please close the previous result CSV file, Output_with_Rxn_no_{ui}.csv, if opened in the background, and rerun again")


get_lambda(gt["chemical_compositions"])

