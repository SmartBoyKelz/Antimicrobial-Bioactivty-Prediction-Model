
import streamlit as st
import pandas as pd
from PIL import Image
import base64
import pickle
import joblib
import os
import subprocess
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

from chembl_webresource_client.new_client import new_client
client = new_client


def calc_Lipinski_desc(smile):
    cmpds = []
    for cmpd in smile:
        cmpdx = Chem.MolFromSmiles(cmpd)

        if cmpdx:
            descriptor = {
                  
        
                  "Molecular wt" : Descriptors.MolWt(cmpdx),
                  "LogP" : Descriptors.MolLogP(cmpdx),
                  "Hydrogen Bond Donors" : Descriptors.NumHDonors(cmpdx),
                  "Hydrogen Bond Donors" : Descriptors.NumHAcceptors(cmpdx),
                  "Aromatic Proportion" : Descriptors.FractionCSP3(cmpdx),
                  "Rotatable bonds" : Descriptors.NumRotatableBonds(cmpdx)
              }

        else:
            descriptor = {
                  
                  "Molecular wt" : "Invalid_cmpd",
                  "LogP" :"Invalid_cmpd",
                  "Hydrogen Bond Donors" : "Invalid_cmpd",
                  "Hydrogen Bond Donors" : "Invalid_cmpd",
                  "Aromatic Proportion" : "Invalid_cmpd",
                  "Rotatable bonds" : "Invalid_cmpd"
              }

        cmpds.append(descriptor)
    return pd.DataFrame(cmpds)


def desc_calc():
    # Performs the descriptor calculation
    bashCommand = "java -Xms2G -Xmx2G -Djava.awt.headless=true -jar ./PaDEL-Descriptor/PaDEL-Descriptor.jar -removesalt -standardizenitro -fingerprints -descriptortypes ./PaDEL-Descriptor/PubchemFingerprinter.xml -dir ./ -file descriptors_output.csv"
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    


def fileDownload(df):
    csv = df.to_csv(index=False)
    b64 = base64.b64encode(csv.encode()).decode()
    href = f'<a href = "data:file/csv;base64, {b64}" download = "prediction.csv">Download Predictions</a>'
    href = f'<a href="data:file/csv;base64,{b64}" download="prediction.csv">Download Predictions</a>'
    return href




def build_model(input_data, load_data):
  #reads in saved regression model
    load_model = joblib.load('XGR.pkl')
    prediction = load_model.predict(input_data)
    st.header('**Predicted Value for Molecule(s) Provided**')
    predicted_output = pd.Series(prediction, name = 'pIC50')
    chmbl_id= pd.Series(load_data[1], name = 'Molecule_chmbl_id')
    molecule_names = [client.molecule.get(chembl_id)['pref_name'] for chembl_id in load_data[1]]
    molecule_names_series = pd.Series(molecule_names, name = 'molecule_name')
    df = pd.concat([molecule_names_series, chmbl_id,df_Mole_descr, predicted_output], axis = 1)
    df.to_csv('pred.smi', sep = '\t', header = True, index = False)
    st.markdown(fileDownload(df), unsafe_allow_html = True)
    st.write("prediction:", prediction)
    st.write("Lipinski's properties:", df)




# Use st.video to display the video
vid = st.video("ui.mp4")
     
#page Title
st.markdown("""

# 🧪 **pIC50 Concentration Against *Staphylococcus aureus* Prediction App** 🧫

This project is dedicated to developing an application that predicts pIC50 concentrations concerning *Staphylococcus aureus*. The pIC50 concentration is a vital metric, calculated as the negative logarithm of the inhibitory concentration. Our app employs advanced predictive models using `XGboost Regression algorithm` to estimate these concentrations, offering valuable insights for drug discovery and therapeutic research focused on combatting *Staphylococcus aureus*. 🌡🔬🦠💊💡

**Credits**
- App built in `Python` + `Streamlit` by  [Aso Kelechi](https://www.linkedin.com/in/seml/) 

## References

1. ChEMBL Database. (Accessed [Date Accessed]) - [https://www.ebi.ac.uk/chembl/](https://www.ebi.ac.uk/chembl/)
2. Data Professor:  [Data Professor](http://youtube.com/dataprofessor)
3. Descriptor calculated using [PaDEL-Descriptor](http://www.yapcwsoft.com/dd/padeldescriptor/) [[Read the Paper]](https://doi.org/10.1002/jcc.21707).
---
""")

with st.sidebar.header('1. Upload CSV data'):
    uploaded_file = st.sidebar.file_uploader("Upload your input file", type=['txt'])
    file_path = "Demo.txt"
    st.sidebar.markdown(f"[Example input file]({file_path})")


if st.sidebar.button('Predict'):
    try:
        load_data = pd.read_table(uploaded_file, sep=" ", header=None)
    except Exception as e:
        st.error(f"An error occurred while reading the file: {str(e)}")
    load_data.to_csv('molecule.smi', sep = '\t', header = False, index = False)
    st.header('**Original input data**')
    st.write(load_data)

    with st.spinner("Calculating descriptors..."):
        desc_calc()
        
    with open('molecule.smi', 'r') as file:
        # Read the file as a CSV with space (' ') as the separator
        dfp = pd.read_csv(file, sep='\t', header=None)
        first_column = dfp.iloc[:, 0]
        second_column = dfp[1]   
    with st.spinner("Calculating Lipinski's Properties..."):
        df_Mole_descr = calc_Lipinski_desc(first_column)
              
    # Read in calculated descriptors and display the dataframe
        st.header('**Calculated molecular descriptors**')
        desc = pd.read_csv('descriptors_output.csv')
        st.write(desc)
        st.write(desc.shape)

    # Read descriptor list used in previously built model
        desc_1 = desc.drop(['Name'], axis = 1)
        desc_subset =desc_1


    # Apply trained model to make prediction on query compounds
        build_model(desc_subset, load_data)
        st.title("3D Molecular Image Visualization")
        smiles_data = load_data[0]
    for i,smiles in enumerate(smiles_data):
            # Generate a molecule object from SMILES
            mol = Chem.MolFromSmiles(smiles)
            molecule_names = [client.molecule.get(chembl_id)['pref_name'] for chembl_id in load_data[1]]
            molecule_name = molecule_names[i]  

            if mol is not None:
                # Render the molecule to an image
                img = Draw.MolToImage(mol, size=(200, 200))
                # Display the image in Streamlit
                st.image(img, caption=f"Molecule: {molecule_name}", use_column_width=False, output_format="PNG")
            else:
                st.write(f"Invalid SMILES: {smiles}")
                
    from PIL import Image
    image = Image.open("computational-biochemicalanalysis-backgrounf.jpg")
    st.image(image, caption='Aso Kelechi 2023', use_column_width=True)
       
else:
    st.info('Upload input data in the sidebar to start!')       
    from PIL import Image
    image = Image.open("computational-biochemicalanalysis-backgrounf.jpg")
    st.image(image, caption='Aso Kelechi 2023', use_column_width=True)

