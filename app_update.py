"""
BioGNN Immunity - Application Streamlit
Interface utilisateur pour la prédiction de propriétés biologiques de molécules
"""

import os
import streamlit as st
import requests
from typing import Tuple, Dict, Optional, Any, List
import streamlit.components.v1 as components

# Import conditionnel de RDKit pour les propriétés moléculaires (optionnel)
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

# Import de PubChemPy
try:
    import pubchempy as pcp
    PUBCHEM_AVAILABLE = True
except ImportError:
    PUBCHEM_AVAILABLE = False

# ============================================================================
# CONFIGURATION
# ============================================================================

st.set_page_config(
    page_title="BioGNN",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# API Configuration
# TODO: INSÉRER L'URL DE VOTRE API GOOGLE CLOUD RUN ICI
GCP_API_URL = "https://your-cloud-run-service.run.app"  # ← REMPLACER PAR L'URL GCP

# Fallback sur secrets.toml si disponible
if 'API_URI' in os.environ:
    BASE_URI = st.secrets.get(os.environ.get('API_URI'), GCP_API_URL)
elif 'cloud_api_uri' in st.secrets:
    BASE_URI = st.secrets['cloud_api_uri']
else:
    BASE_URI = GCP_API_URL

# Assurer que l'URL se termine par '/'
BASE_URI = BASE_URI if BASE_URI.endswith('/') else BASE_URI + '/'

# ============================================================================
# STYLES CSS PERSONNALISÉS
# ============================================================================

st.markdown("""
<style>
    /* Style principal */
    .main {
        background-color: #2d3a2e;
    }

    /* Cartes avec style olive/vert */
    .prediction-card {
        background: linear-gradient(135deg, #4a5d4e 0%, #5a6d5e 100%);
        padding: 2rem;
        border-radius: 15px;
        box-shadow: 0 4px 6px rgba(0, 0, 0, 0.3);
        margin: 1rem 0;
        border: 2px solid #6b7d6b;
    }

    .result-card {
        background: linear-gradient(135deg, #5a6d5e 0%, #6a7d6e 100%);
        padding: 1.5rem;
        border-radius: 12px;
        margin: 0.5rem 0;
        border-left: 4px solid #9db89d;
    }

    /* Texte de prédiction */
    .prediction-text {
        color: #b8e986;
        font-size: 2rem;
        font-weight: 700;
        text-align: center;
        margin: 1rem 0;
        text-shadow: 2px 2px 4px rgba(0, 0, 0, 0.3);
    }

    .prediction-label {
        color: #d4d4d4;
        font-size: 0.9rem;
        font-weight: 500;
        text-transform: uppercase;
        letter-spacing: 1px;
        opacity: 0.8;
    }

    /* Input boxes */
    .stTextInput > div > div > input {
        background-color: #3d4a3e;
        color: #e0e0e0;
        border: 2px solid #5a6d5e;
        border-radius: 8px;
    }

    .stSelectbox > div > div > select {
        background-color: #3d4a3e;
        color: #e0e0e0;
        border: 2px solid #5a6d5e;
    }

    /* Boutons */
    .stButton > button {
        background: linear-gradient(135deg, #6b8e6b 0%, #8aae8a 100%);
        color: white;
        font-weight: 600;
        border: none;
        border-radius: 8px;
        padding: 0.75rem 2rem;
        font-size: 1.1rem;
        box-shadow: 0 4px 6px rgba(0, 0, 0, 0.2);
        transition: all 0.3s ease;
    }

    .stButton > button:hover {
        background: linear-gradient(135deg, #7a9d7a 0%, #9bbd9b 100%);
        box-shadow: 0 6px 8px rgba(0, 0, 0, 0.3);
        transform: translateY(-2px);
    }

    /* Titre principal */
    h1 {
        color: #b8e986;
        text-shadow: 2px 2px 4px rgba(0, 0, 0, 0.3);
    }

    h2, h3 {
        color: #d4d4d4;
    }

    /* Molécule display */
    .molecule-container {
        background-color: #f8f8f8;
        padding: 1rem;
        border-radius: 10px;
        margin: 1rem 0;
        box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
    }

    /* Info boxes */
    .info-box {
        background-color: #4a5d4e;
        padding: 1rem;
        border-radius: 8px;
        border-left: 4px solid #b8e986;
        color: #e0e0e0;
        margin: 0.5rem 0;
    }

    /* Sidebar */
    .css-1d391kg {
        background-color: #252f26;
    }

    /* Expander */
    .streamlit-expanderHeader {
        background-color: #3d4a3e;
        border-radius: 8px;
        color: #e0e0e0;
    }

    /* Métrics */
    .metric-container {
        background: linear-gradient(135deg, #4a5d4e 0%, #5a6d5e 100%);
        padding: 1rem;
        border-radius: 10px;
        text-align: center;
        margin: 0.5rem 0;
    }
</style>
""", unsafe_allow_html=True)

# ============================================================================
# FONCTIONS UTILITAIRES
# ============================================================================

def validate_smiles(smiles: str) -> Tuple[bool, str]:
    """
    Valide un SMILES et retourne (is_valid, message)
    """
    if not RDKIT_AVAILABLE:
        # Validation basique sans RDKit
        if not smiles or len(smiles) < 1:
            return False, "SMILES vide"
        # Vérification basique de caractères valides
        valid_chars = set("CNOPSFClBrI[]()=#@+-0123456789cnops")
        if not all(c in valid_chars for c in smiles):
            return False, "Caractères invalides dans le SMILES"
        return True, "SMILES accepté (validation basique) ✓"

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "SMILES invalide - structure non reconnue"
        return True, "SMILES valide ✓"
    except Exception as e:
        return False, f"Erreur: {str(e)}"

def render_molecule_3d(smiles: str, height: int = 400, width: int = 400) -> str:
    """
    Génère une visualisation 3D interactive d'une molécule avec py3Dmol

    Args:
        smiles: String SMILES de la molécule
        height: Hauteur du viewer en pixels
        width: Largeur du viewer en pixels

    Returns:
        HTML string pour affichage avec st.components.html
    """
    html = f"""
    <html>
    <head>
        <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
    </head>
    <body>
        <div id="container" style="height: {height}px; width: {width}px; position: relative;"></div>
        <script>
            let element = document.getElementById('container');
            let config = {{ backgroundColor: 'white' }};
            let viewer = $3Dmol.createViewer(element, config);

            // Convertir SMILES en structure 3D via PubChem
            // Fallback: afficher juste le SMILES si la conversion échoue
            fetch('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/SDF')
                .then(response => response.text())
                .then(sdf => {{
                    viewer.addModel(sdf, 'sdf');
                    viewer.setStyle({{}}, {{stick: {{}}, sphere: {{radius: 0.3}}}});
                    viewer.zoomTo();
                    viewer.render();
                }})
                .catch(error => {{
                    // Si PubChem échoue, afficher un message
                    element.innerHTML = '<div style="padding: 20px; text-align: center; color: #666;">' +
                        '<p><strong>Molécule:</strong> {smiles}</p>' +
                        '<p style="font-size: 0.9em;">Visualisation 3D non disponible</p>' +
                        '</div>';
                }});
        </script>
    </body>
    </html>
    """
    return html

def get_pubchem_publications(smiles: str, max_results: int = 5) -> List[Dict[str, Any]]:
    """
    Récupère les publications scientifiques liées à une molécule via PubChem

    Args:
        smiles: String SMILES de la molécule
        max_results: Nombre maximum de publications à retourner

    Returns:
        Liste de dictionnaires contenant les informations des publications
    """
    if not PUBCHEM_AVAILABLE:
        return [{
            "error": "PubChemPy non disponible",
            "message": "Impossible de récupérer les publications"
        }]

    try:
        # Rechercher le composé via SMILES
        compounds = pcp.get_compounds(smiles, 'smiles')

        if not compounds:
            return [{
                "title": "Aucun composé trouvé",
                "message": "Cette molécule n'est pas référencée dans PubChem"
            }]

        compound = compounds[0]
        cid = compound.cid

        # Récupérer les publications via l'API PubChem REST
        # Utiliser l'API PubMed via PubChem
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/xrefs/PubMedID/JSON"

        response = requests.get(url, timeout=10)

        if response.status_code != 200:
            return [{
                "title": f"Composé trouvé: {compound.iupac_name or 'Sans nom'}",
                "cid": cid,
                "message": "Aucune publication PubMed trouvée pour ce composé"
            }]

        data = response.json()
        pmids = data.get('InformationList', {}).get('Information', [{}])[0].get('PubMedID', [])

        if not pmids:
            return [{
                "title": f"Composé: {compound.iupac_name or 'CID ' + str(cid)}",
                "cid": cid,
                "message": "Aucune publication disponible"
            }]

        # Limiter au nombre demandé
        pmids = pmids[:max_results]

        # Récupérer les détails des publications via PubMed
        publications = []
        for pmid in pmids:
            try:
                # API PubMed pour récupérer les détails
                pubmed_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id={pmid}&retmode=json"
                pub_response = requests.get(pubmed_url, timeout=5)

                if pub_response.status_code == 200:
                    pub_data = pub_response.json()
                    result = pub_data.get('result', {}).get(str(pmid), {})

                    publications.append({
                        "pmid": pmid,
                        "title": result.get('title', 'Titre non disponible'),
                        "authors": ', '.join([author.get('name', '') for author in result.get('authors', [])[:3]]),
                        "journal": result.get('fulljournalname', result.get('source', 'Journal inconnu')),
                        "year": result.get('pubdate', 'Date inconnue').split()[0],
                        "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
                    })
            except Exception:
                continue

        return publications if publications else [{
            "title": "Erreur lors de la récupération des publications",
            "message": "Impossible de charger les détails depuis PubMed"
        }]

    except Exception as e:
        return [{
            "error": "Erreur",
            "message": f"Impossible de récupérer les publications: {str(e)}"
        }]

def get_organism_info(organism_name: str) -> str:
    """
    Retourne les informations sur un organisme modèle à partir d'un dictionnaire

    Args:
        organism_name: Nom de l'organisme à rechercher

    Returns:
        Description de l'organisme
    """
    organism_info = {
        'Rattus norvegicus': "Le rat norvégien (Rattus norvegicus), également appelé rat brun ou surmulot, est une espèce de rongeurs très utilisée en recherche biomédicale. C'est un organisme modèle important pour l'étude de la physiologie, de la pharmacologie et du comportement. Sa proximité génétique avec l'homme en fait un modèle précieux pour la recherche sur les maladies humaines.",

        'Saccharomyces cerevisiae': "Saccharomyces cerevisiae, communément appelée levure de boulanger ou levure de bière, est un organisme eucaryote unicellulaire. C'est l'un des organismes modèles les plus étudiés en biologie cellulaire et moléculaire. Son génome a été le premier d'un eucaryote à être entièrement séquencé en 1996.",

        'Bos taurus': "Bos taurus est le nom scientifique du bœuf domestique ou bovin. Bien que principalement élevé pour l'agriculture, il est également utilisé en recherche biomédicale, notamment pour l'étude du métabolisme, de la reproduction et de certaines maladies zoonotiques.",

        'Cavia porcellus': "Le cobaye (Cavia porcellus), également appelé cochon d'Inde, est un rongeur originaire d'Amérique du Sud. C'est un organisme modèle important en immunologie, en nutrition et en recherche sur les maladies infectieuses. Contrairement à la plupart des mammifères, il ne peut pas synthétiser la vitamine C.",

        'Oryctolagus cuniculus': "Le lapin de garenne (Oryctolagus cuniculus) est utilisé comme organisme modèle en recherche biomédicale, notamment en immunologie, en toxicologie et en recherche cardiovasculaire. Il est particulièrement utile pour la production d'anticorps et l'étude de certaines maladies oculaires.",

        'Cricetulus griseus': "Le hamster chinois (Cricetulus griseus) est un petit rongeur originaire de Chine et de Mongolie. Ses cellules ovariennes (cellules CHO) sont largement utilisées en biotechnologie pour la production de protéines thérapeutiques et en recherche sur la génétique cellulaire.",

        'Chlorocebus sabaeus': "Le singe vert (Chlorocebus sabaeus) est un primate originaire d'Afrique de l'Ouest. Il est utilisé comme organisme modèle en recherche biomédicale, notamment pour l'étude des maladies infectieuses, du VIH/SIDA et des troubles neurologiques en raison de sa proximité phylogénétique avec l'homme.",

        'Mus musculus': "La souris de laboratoire (Mus musculus) est l'organisme modèle mammifère le plus utilisé en recherche biomédicale. Son génome, très similaire à celui de l'homme, permet d'étudier de nombreuses maladies humaines. Elle est utilisée dans pratiquement tous les domaines de la recherche biologique.",

        'Escherichia coli': "Escherichia coli (E. coli) est une bactérie intestinale commune qui constitue l'un des organismes modèles les plus importants en biologie moléculaire et en génétique microbienne. Elle est largement utilisée en biotechnologie pour la production de protéines recombinantes et comme hôte pour le clonage génétique.",

        'Trypanosoma cruzi': "Trypanosoma cruzi est un protozoaire parasite responsable de la maladie de Chagas, une maladie tropicale négligée affectant des millions de personnes en Amérique latine. C'est un organisme modèle important pour l'étude des maladies parasitaires et le développement de nouveaux traitements antiparasitaires.",

        'Human immunodeficiency virus 1': "Le virus de l'immunodéficience humaine de type 1 (VIH-1) est le rétrovirus responsable du SIDA. Il est intensivement étudié en virologie et en immunologie pour comprendre les mécanismes d'infection virale et développer des thérapies antivirales et des vaccins.",

        'Mycobacterium tuberculosis': "Mycobacterium tuberculosis est la bactérie responsable de la tuberculose, l'une des maladies infectieuses les plus mortelles au monde. C'est un organisme modèle crucial pour la recherche sur les maladies infectieuses, le développement d'antibiotiques et l'étude des mécanismes de résistance aux médicaments.",

        'Drosophila': "La drosophile (Drosophila melanogaster), ou mouche du vinaigre, est un organisme modèle fondamental en génétique et en biologie du développement. Son cycle de vie court, sa facilité de culture et son génome bien caractérisé en font un outil précieux pour étudier l'hérédité, le développement embryonnaire et les maladies neurodégénératives.",

        'Schistosoma mansoni': "Schistosoma mansoni est un ver parasite responsable de la schistosomiase, une maladie tropicale affectant des millions de personnes. C'est un organisme modèle important pour l'étude des helminthiases et le développement de médicaments antiparasitaires.",

        'Caenorhabditis elegans': "Caenorhabditis elegans est un nématode transparent d'environ 1 mm de long. C'est un organisme modèle majeur en biologie du développement, en génétique et en neurobiologie. Son système nerveux simple mais complet, composé de seulement 302 neurones, est entièrement cartographié.",

        'Bacillus subtilis': "Bacillus subtilis est une bactérie Gram-positive couramment trouvée dans le sol. C'est un organisme modèle important en microbiologie pour l'étude de la formation des spores, de la différenciation cellulaire et de la production d'enzymes industrielles.",

        'Leishmania mexicana': "Leishmania mexicana est un protozoaire parasite responsable de la leishmaniose cutanée. C'est un organisme modèle pour l'étude des maladies parasitaires tropicales et le développement de nouveaux traitements contre les leishmanioses.",

        'Trypanosoma brucei TREU927': "Trypanosoma brucei est un protozoaire parasite responsable de la maladie du sommeil en Afrique. La souche TREU927 est utilisée comme organisme modèle de référence pour l'étude de la biologie des trypanosomes et le développement de médicaments antiparasitaires."
    }

    # Retourner l'information ou un message par défaut
    return organism_info.get(organism_name, f"Organisme modèle utilisé en recherche biomédicale pour l'étude de diverses propriétés biologiques et le développement de nouveaux composés thérapeutiques.")

def call_api(smiles: str, property_name: str, organism: str = "Homo sapiens", model: str = "immunity") -> Dict[str, Any]:
    """
    Appelle l'API GCP pour obtenir une prédiction

    Args:
        smiles: SMILES de la molécule
        property_name: Propriété biologique à prédire
        organism: Organisme cible pour la prédiction
        model: Modèle à utiliser (immunity ou antiox)

    TODO: Adapter les paramètres selon votre API
    """
    try:
        url = BASE_URI + "predict"

        # Paramètres de la requête - ADAPTER SELON VOTRE API
        params = {
            "smiles": smiles,
            "property": property_name,
            "organism": organism,
            "model": model
        }

        # Appel API
        response = requests.get(url, params=params, timeout=30)

        if response.status_code == 200:
            return {
                "success": True,
                "data": response.json(),
                "status_code": 200
            }
        else:
            return {
                "success": False,
                "error": f"Erreur API: Status {response.status_code}",
                "status_code": response.status_code
            }

    except requests.exceptions.Timeout:
        return {
            "success": False,
            "error": "Timeout - L'API ne répond pas",
            "status_code": 408
        }
    except requests.exceptions.ConnectionError:
        return {
            "success": False,
            "error": "Impossible de se connecter à l'API",
            "status_code": 503
        }
    except Exception as e:
        return {
            "success": False,
            "error": f"Erreur inattendue: {str(e)}",
            "status_code": 500
        }

# ============================================================================
# INTERFACE PRINCIPALE
# ============================================================================

def main():
    # En-tête
    st.markdown("<h1 style='text-align: center;'>🧬 BioGNN 🧬 <br>From atoms to action</h1>",unsafe_allow_html=True)
    st.markdown("<p style='text-align: center; color: #b0b0b0; font-size: 1.2rem;'>Prédiction de propriétés biologiques par Graph Neural Networks</p>", unsafe_allow_html=True)

    st.markdown("---")

    # Sidebar - Informations et paramètres
    with st.sidebar:
        # Exemples de SMILES
        st.markdown("### 📋 Exemples de SMILES")
        examples = {
            "Aspirine": "CC(=O)Oc1ccccc1C(=O)O",
            "Caféine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "Éthanol": "CCO",
            "Glucose": "C(C1C(C(C(C(O1)O)O)O)O)O",
            "Pénicilline G": "CC1(C(N2C(S1)C(C2=O)NC(=O)Cc3ccccc3)C(=O)O)C"
        }

        for name, smiles in examples.items():
            if st.button(f"💊 {name}", use_container_width=True):
                st.session_state['example_smiles'] = smiles

        st.markdown("---")
        st.markdown("### 📖 À propos")
        st.info("""
        **BioGNN** utilise des Graph Neural Networks pour prédire
        les propriétés biologiques de molécules à partir de leur structure SMILES et d'un organisme modèle cible.
        """)

    # Zone principale - Input
    col1, col2 = st.columns([1, 1])

    with col1:
        st.markdown("#### 🧪 Molécule (SMILES)")

        # Input SMILES
        default_smiles = st.session_state.get('example_smiles', '')
        smiles_input = st.text_input(
            "Entrez le SMILES de votre molécule:",
            value=default_smiles,
            placeholder="CC(=O)Oc1ccccc1C(=O)O",
            label_visibility="collapsed"
        )

        # Validation
        if smiles_input:
            is_valid, msg = validate_smiles(smiles_input)
            if is_valid:
                st.success(msg)
            else:
                st.error(msg)

    with col2:
        st.markdown("#### 🐁 Organisme modèle ")

        # Sélection du modèle
        selected_organism = st.selectbox(
            "Sélectionnez le modèle:",
            [ 'Rattus norvegicus','Saccharomyces cerevisiae','Bos taurus','Cavia porcellus','Oryctolagus cuniculus','Cricetulus griseus','Chlorocebus sabaeus','Mus musculus','Escherichia coli','Trypanosoma cruzi','Human immunodeficiency virus 1','Mycobacterium tuberculosis','Drosophila','Schistosoma mansoni','Caenorhabditis elegans','Bacillus subtilis','Leishmania mexicana','Trypanosoma brucei TREU927'],
            label_visibility="collapsed",
            index=7
        )

    # Bouton de prédiction
    st.markdown(
    """
    <style>
    button[kind="primary"] {
        height: 70px;
        font-weight: bold;
    }

    button[kind="primary"] * {
        font-size: 25px;
    }
    </style>
    """,
    unsafe_allow_html=True)

    st.markdown("<br>", unsafe_allow_html=True)
    col_btn1, col_btn2, col_btn3 = st.columns([1, 2, 1])
    with col_btn2:
        predict_button = st.button("PRÉDIRE", use_container_width=True,type="primary")

    st.markdown("---")

    # ========================================================================
    # AFFICHAGE DES RÉSULTATS
    # ========================================================================

    if predict_button and smiles_input:
        is_valid, _ = validate_smiles(smiles_input)

        if not is_valid:
            st.error("⚠️ Veuillez entrer un SMILES valide avant de prédire")
        else:
            # Affichage de la molécule
            st.markdown("<h3 style='text-align: center;'>La molécule</h3>",unsafe_allow_html=True)

            col_mol1, col_mol2 = st.columns([1, 1])

            with col_mol1:
                st.markdown("#### 🧬 Structure 3D Interactive")
                # Générer la visualisation 3D avec py3Dmol
                mol_html = render_molecule_3d(smiles_input, height=400, width=500)
                components.html(mol_html, height=450, scrolling=False)

            with col_mol2:
                st.markdown("#### 📚 Publications Scientifiques")

                with st.spinner("🔍 Recherche de publications..."):
                    publications = get_pubchem_publications(smiles_input, max_results=2)

                if publications and not publications[0].get('error'):
                    for i, pub in enumerate(publications, 1):
                        if 'pmid' in pub:
                            # Publication complète
                            st.markdown(f"""
                            <div class="info-box" style="margin-bottom: 1rem;">
                                <p style="margin: 0; font-size: 0.85rem; color: #9db89d; font-weight: 600;">
                                    #{i} • {pub['year']}
                                </p>
                                <p style="margin: 0.5rem 0; font-weight: 600; font-size: 0.95rem;">
                                    {pub['title']}
                                </p>
                                <p style="margin: 0.3rem 0; font-size: 0.85rem; color: #c0c0c0;">
                                    <strong>Auteurs:</strong> {pub['authors']}
                                </p>
                                <p style="margin: 0.3rem 0; font-size: 0.85rem; color: #c0c0c0;">
                                    <strong>Journal:</strong> {pub['journal']}
                                </p>
                                <a href="{pub['url']}" target="_blank" style="
                                    display: inline-block;
                                    margin-top: 0.5rem;
                                    padding: 0.3rem 0.8rem;
                                    background: #6b8e6b;
                                    color: white;
                                    text-decoration: none;
                                    border-radius: 5px;
                                    font-size: 0.85rem;
                                    transition: background 0.3s;
                                ">
                                    📖 Lire sur PubMed
                                </a>
                            </div>
                            """, unsafe_allow_html=True)
                        else:
                            # Message d'information ou d'erreur
                            st.info(f"ℹ️ {pub.get('title', 'Information')} - {pub.get('message', '')}")
                else:
                    st.warning("⚠️ Aucune publication trouvée pour cette molécule dans PubChem/PubMed")

            st.markdown("---")

            # Affichage de l'espèce modèle
            st.markdown("<h3 style='text-align: center;'>L'espèce modèle</h3>",unsafe_allow_html=True)

            col_species1, col_species2, col_species3 = st.columns([1, 2, 1])

            with col_species2:
                # Récupérer les informations sur l'organisme
                organism_description = get_organism_info(selected_organism)

                st.markdown(f"""
                <div class="prediction-card" style="text-align: center;">
                    <p class="prediction-label">ORGANISME SÉLECTIONNÉ</p>
                    <p class="prediction-text" style="font-size: 1.8rem;">{selected_organism}</p>
                    <div class="info-box" style="margin-top: 1rem; text-align: left;">
                        <p style="margin: 0.5rem 0; color: #e0e0e0;">
                            {organism_description}
                        </p>
                    </div>
                </div>
                """, unsafe_allow_html=True)

            st.markdown("---")

            # Appel à l'API
            st.markdown("<h3 style='text-align: center;'>Prédiction</h3>",unsafe_allow_html=True)

            # Afficher les paramètres de prédiction
            st.info(f"🤖 **Modèle:** {model_choice} | 🧬 **Organisme:** {selected_organism} | 🎯 **Propriété:** {selected_property}")

            with st.spinner("🔄 Analyse en cours..."):
                result = call_api(smiles_input, selected_property, selected_organism, model_choice)

            if result["success"]:
                # TODO: ADAPTER L'AFFICHAGE SELON LA STRUCTURE DE VOTRE API

                st.markdown('<div class="prediction-card">', unsafe_allow_html=True)

                st.markdown('<p class="prediction-label">PRÉDICTION</p>', unsafe_allow_html=True)
                st.markdown('<p class="prediction-text">Candidat prometteur</p>', unsafe_allow_html=True)

                # Afficher les données brutes (à adapter)
                with st.expander("📊 Détails de la prédiction"):
                    st.json(result["data"])

                st.markdown('</div>', unsafe_allow_html=True)

                # Métriques supplémentaires (à adapter selon votre API)
                st.markdown("#### 📈 Métriques")
                metric_cols = st.columns(3)

                with metric_cols[0]:
                    st.markdown('<div class="metric-container">', unsafe_allow_html=True)
                    st.metric("Score de confiance", "85%")  # TODO: Remplacer par vraie valeur
                    st.markdown('</div>', unsafe_allow_html=True)

                with metric_cols[1]:
                    st.markdown('<div class="metric-container">', unsafe_allow_html=True)
                    st.metric("Probabilité", "0.78")  # TODO: Remplacer par vraie valeur
                    st.markdown('</div>', unsafe_allow_html=True)

                with metric_cols[2]:
                    st.markdown('<div class="metric-container">', unsafe_allow_html=True)
                    st.metric("Classe prédite", "Actif")  # TODO: Remplacer par vraie valeur
                    st.markdown('</div>', unsafe_allow_html=True)

            else:
                st.error(f"❌ {result['error']}")
                st.info("💡 Vérifiez que l'URL de l'API est correctement configurée dans le code")

    elif predict_button:
        st.warning("⚠️ Veuillez entrer un SMILES valide")

    # Footer
    st.markdown("---")
    st.markdown("""
    <div style='text-align: center; color: #808080; padding: 2rem;'>
        <p>🧬 BioGNN - Propulsé par Graph Neural Networks</p>
    </div>
    """, unsafe_allow_html=True)

# ============================================================================
# POINT D'ENTRÉE
# ============================================================================

if __name__ == "__main__":
    main()
