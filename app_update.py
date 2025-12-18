"""
BioGNN Immunity - Application Streamlit
Interface utilisateur pour la prédiction de propriétés biologiques de molécules
"""

import os
import streamlit as st
import requests
from typing import Tuple, Dict, Optional, Any, List
import streamlit.components.v1 as components

# ===============================
# Mapping des organismes supportés pour la toxicité
# ===============================
ORGANISMS_TOXICITY_MAPPING = {
    'rattus': 'Rattus norvegicus',
    'equus': 'Equus caballus',
    'h1n1': 'Influenza A virus (H1N1)',
    'ecoli': 'Escherichia coli',
    'hiv': 'Human immunodeficiency virus 1',
    'gondii': 'Toxoplasma gondii'
}

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
BASE_URI = "https://biognn-third-api-223608804946.europe-west1.run.app/"

# def get_api_url():
#     """
#     Détermine l'URL de l'API à utiliser selon l'environnement:
#     1. Si une variable d'environnement API_URL est définie, l'utiliser
#     2. Si en local (détection via secrets), utiliser local_api_uri ou local_docker_uri
#     3. Si cloud_api_uri est défini dans secrets, l'utiliser
#     4. Sinon, utiliser l'URL GCP par défaut
#     """
#     # Priority 1: Variable d'environnement directe
#     if 'API_URL' in os.environ:
#         url = os.environ['API_URL']
#         return url if url.endswith('/') else url + '/'

#     # Priority 2: Vérifier si secrets.toml existe et contient des URLs
#     try:
#         # Détecter si on est en environnement local
#         # Si local_api_uri ou local_docker_uri existe dans secrets, on est probablement en local
#         if hasattr(st.secrets, 'get'):
#             # Mode local: priorité à local_api_uri si disponible
#             if 'local_api_uri' in st.secrets:
#                 # Vérifier si on force l'usage de l'API locale
#                 use_local = st.secrets.get('use_local_api', False)
#                 if use_local:
#                     url = st.secrets['local_api_uri']
#                     return url if url.endswith('/') else url + '/'

#             # Docker local
#             if 'local_docker_uri' in st.secrets:
#                 use_docker = st.secrets.get('use_docker_api', False)
#                 if use_docker:
#                     url = st.secrets['local_docker_uri']
#                     return url if url.endswith('/') else url + '/'

#             # Cloud API (production ou ancien déploiement)
#             if 'cloud_api_uri' in st.secrets:
#                 url = st.secrets['cloud_api_uri']
#                 return url if url.endswith('/') else url + '/'
#     except Exception:
#         # Si erreur d'accès aux secrets, continuer avec le fallback
#         pass

#     # Priority 3: URL GCP par défaut
#     return GCP_API_URL if GCP_API_URL.endswith('/') else GCP_API_URL + '/'

# BASE_URI = get_api_url()

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

    /* Team Modal */
    .team-modal {
        position: fixed;
        z-index: 9999;
        left: 0;
        top: 0;
        width: 100%;
        height: 100%;
        overflow: auto;
        background-color: rgba(0, 0, 0, 0.8);
        display: flex;
        align-items: center;
        justify-content: center;
    }

    .team-modal-content {
        background: linear-gradient(135deg, #2d3a2e 0%, #3d4a3e 100%);
        padding: 2rem;
        border-radius: 15px;
        max-width: 95%;
        width: 1200px;
        max-height: 90vh;
        overflow-y: auto;
        box-shadow: 0 8px 16px rgba(0, 0, 0, 0.5);
        border: 2px solid #6b8e6b;
    }

    .team-grid {
        display: flex;
        flex-direction: row;
        gap: 1.5rem;
        margin-top: 1.5rem;
        justify-content: center;
        flex-wrap: wrap;
    }

    .team-member {
        background: linear-gradient(135deg, #4a5d4e 0%, #5a6d5e 100%);
        padding: 1.5rem;
        border-radius: 12px;
        border: 2px solid #6b7d6b;
        text-align: center;
        flex: 0 1 250px;
        min-width: 220px;
    }

    .team-member img {
        width: 150px;
        height: 150px;
        min-width: 150px;
        min-height: 150px;
        max-width: 150px;
        max-height: 150px;
        border-radius: 50%;
        object-fit: cover;
        object-position: center;
        border: 3px solid #b8e986;
        margin-bottom: 1rem;
        display: block;
        margin-left: auto;
        margin-right: auto;
    }

    .team-member h3 {
        color: #b8e986;
        margin: 0.5rem 0;
        font-size: 1.1rem;
        line-height: 1.3;
    }

    .team-member p {
        color: #d4d4d4;
        font-size: 0.85rem;
        margin: 0.5rem 0;
        line-height: 1.4;
    }

    .team-links {
        display: flex;
        justify-content: center;
        gap: 0.8rem;
        margin-top: 0.8rem;
    }

    .team-links a {
        display: inline-block;
        padding: 0.4rem 0.8rem;
        background: #6b8e6b;
        color: white;
        text-decoration: none;
        border-radius: 5px;
        font-size: 0.85rem;
        transition: background 0.3s;
    }

    .team-links a:hover {
        background: #7a9d7a;
    }

    .close-modal {
        color: #b8e986;
        float: right;
        font-size: 2rem;
        font-weight: bold;
        cursor: pointer;
        line-height: 1;
    }

    .close-modal:hover {
        color: #9db89d;
    }
</style>
""", unsafe_allow_html=True)

# ============================================================================
# FONCTIONS UTILITAIRES
# ============================================================================

@st.dialog("🧬 Notre Équipe", width="large")
def show_team_modal():
    """
    Affiche un modal avec les informations des membres de l'équipe
    """
    team_members = [
        {
            "name": "Victor Carré",
            "title": "PhD in Organic Chemistry\nData Scientist",
            "photo": "https://media.licdn.com/dms/image/v2/D4E03AQGVWlViiqc8YA/profile-displayphoto-shrink_800_800/profile-displayphoto-shrink_800_800/0/1721578462243?e=1767225600&v=beta&t=PXmCdebrZyzU2R3SLz_0VEEkzK-2uOXEn8pLOaXJf_M",
            "description": "Bringing photochemistry expertise with Data Science to decode molecular behavior.",
            "linkedin": "https://www.linkedin.com/in/victor-carré",
            "github": "https://github.com/victorcarre6"
        },
        {
            "name": "Nisha Dwivedi",
            "title": "PhD in Bioinformatics\nData Scientist",
            "photo": "https://avatars.githubusercontent.com/u/97964928?v=4",
            "description": "From decoding genomes to decoding data,\nI bridge biology and data science.",
            "linkedin": "https://www.linkedin.com/in/nisha-dwivedi-108b64206/",
            "github": "https://github.com/nishadwivedi97"
        },
        {
            "name": "Jalil Kheloufi",
            "title": "Data Scientist\n ",
            "photo": "https://media.licdn.com/dms/image/v2/D4E03AQEzs0-wkrE4gg/profile-displayphoto-shrink_800_800/B4EZQ8udQ9G4Ac-/0/1736185599605?e=1767225600&v=beta&t=7oR82b7G8SmXhePojlejLQXaqzdV4n1VmfELTReKOzk",
            "description": "CRM Developer at Salesforce.\nPassionate about data science and AI.",
            "linkedin": "https://www.linkedin.com/in/jalilkheloufi/",
            "github": "https://github.com/Soipadeg"
        },
        {
            "name": "Jean-Charles Bodart",
            "title": "Data Scientist\n ",
            "photo": "https://media.licdn.com/dms/image/v2/C4E03AQHPxvytYnRNVQ/profile-displayphoto-shrink_800_800/profile-displayphoto-shrink_800_800/0/1659509572655?e=1767225600&v=beta&t=fem0yJcUaLU4O4CLm8sp8Wh9yYzhgJjQI6-53rXOG5g",
            "description": "Turning raw data into strategic decisions through data-driven thinking.",
            "linkedin": "https://www.linkedin.com/in/jean-charles-bodart-492a40a0/",
            "github": "https://github.com/jeancharlesbodart-commits"
        }
    ]

    # Afficher les membres en colonnes
    cols = st.columns(4)

    for idx, member in enumerate(team_members):
        with cols[idx]:
            # Photo
            st.markdown(f"""
                <div style="text-align: center;">
                    <img src="{member['photo']}"
                         style="width: 150px; height: 150px; border-radius: 50%;
                                object-fit: cover; border: 3px solid #b8e986;
                                margin-bottom: 1rem;">
                </div>
            """, unsafe_allow_html=True)

            # Nom et titre avec hauteur fixe
            st.markdown(f"<h4 style='text-align: center; color: #b8e986; font-size: 1.5rem; margin: 0.5rem 0;'>{member['name']}</h4>", unsafe_allow_html=True)
            st.markdown(f"<div style='text-align: center; color: #9db89d; font-size: 1rem; margin: 0.3rem 0; height: 2.5rem; line-height: 1.25rem; white-space: pre-line;'>{member['title']}</div>", unsafe_allow_html=True)

            # Description avec hauteur fixe
            st.markdown(f"<div style='text-align: center; color: #d4d4d4; font-size: 1rem; margin: 1rem 0; height: 3rem;'>{member['description']}</div>", unsafe_allow_html=True)

            # Liens
            col1, col2 = st.columns(2)
            with col1:
                st.markdown(f"<a href='{member['linkedin']}' target='_blank' style='display: block; text-align: center; padding: 0.5rem; background: #0077b5; color: white; text-decoration: none; border-radius: 5px; font-size: 0.85rem;'>LinkedIn</a>", unsafe_allow_html=True)
            with col2:
                st.markdown(f"<a href='{member['github']}' target='_blank' style='display: block; text-align: center; padding: 0.5rem; background: #333; color: white; text-decoration: none; border-radius: 5px; font-size: 0.85rem;'>GitHub</a>", unsafe_allow_html=True)

@st.dialog("🧠 Informations sur le Modèle", width="large")
def show_model_info():
    """
    Affiche un modal avec les informations techniques sur le modèle BioGNN
    """
    model_info = [
        {
            "title": "Source des données",
            "logo": "https://avatars.githubusercontent.com/u/3062531?s=280&v=4",
            "description": "Les données d'entrainement proviennent de la base publique ChemBL, regroupant plus de 2.8 millions de tests biologiques et 1.7 millions de molécules.",
            "ref": "https://www.ebi.ac.uk/chembl/"
        },
        {
            "title": "Neural Networks",
            "logo": "https://thumbs.dreamstime.com/b/neural-network-cloud-technologies-global-database-artificial-intelligence-bright-black-white-background-big-data-d-200490439.jpg",
            "description": "BioGNN est un modèle de Deep Learning reposant sur une architecture de réseaux de neurones graphiques et des mécanismes d'attention head inspiré de l'IA générative.",
            "ref": "https://pytorch.org/docs/stable/nn.html"
        },
        {
            "title": "Technologies",
            "logo": "https://avatars.githubusercontent.com/u/57251745?s=280&v=4",
            "description": "Le projet se base sur le paquet PyTorch, exploite la librairie de chemoinformatique RDKit, et le framework Optuna pour l'optimisation des réseaux neuronaux.",
            "ref": "https://www.rdkit.org/",
            "ref2": "https://optuna.org/"
        },
        {
            "title": "Performances",
            "logo": "https://cdn.dribbble.com/userupload/41870104/file/original-2b633babbb87748069a3c2924124a344.gif",
            "description": "Le modèle atteint une précision de 85% sur un jeu de test indépendant, montrant sa capacité à généraliser sur des molécules. Les limitations du modèle se situent dans sa capacité à généraliser sur le type d'espèce.",
            "ref" : "https://arxiv.org/abs/2507.03430"
        }
    ]

    # Afficher les cartes en colonnes (4 colonnes maintenant)
    cols = st.columns(4)

    for idx, info in enumerate(model_info):
        with cols[idx]:
            # Logo avec hauteur fixe
            if 'logo' in info and info['logo']:
                st.markdown(f"""
                    <div style="text-align: center; height: 140px; display: flex; align-items: center; justify-content: center;">
                        <img src="{info['logo']}"
                             style="width: 120px; height: 120px; border-radius: 10px;
                                    object-fit: cover; border: 2px solid #b8e986;">
                    </div>
                """, unsafe_allow_html=True)
            else:
                st.markdown("<div style='height: 140px;'></div>", unsafe_allow_html=True)

            # Titre avec hauteur fixe
            st.markdown(f"<div style='text-align: center; color: #b8e986; font-size: 1.3rem; margin: 0.5rem 0; height: 2.6rem; display: flex; align-items: center; justify-content: center; font-weight: bold;'>{info['title']}</div>", unsafe_allow_html=True)

            # Description avec hauteur fixe
            st.markdown(f"<div style='text-align: center; color: #d4d4d4; font-size: 1.1rem; margin: 1rem 0; line-height: 1.4; height: 100px; display: flex; align-items: top; justify-content: center;'>{info['description']}</div>", unsafe_allow_html=True)

            # Liens (gérer ref et ref2 pour Librairies) avec hauteur fixe
            st.markdown("<div style='height: 50px; display: flex; align-items: center;'>", unsafe_allow_html=True)
            if 'ref2' in info:
                # Cas spécial: deux liens (Librairies)
                col1, col2 = st.columns(2)
                with col1:
                    st.markdown(f"<a href='{info['ref']}' target='_blank' style='display: block; text-align: center; padding: 0.5rem; background: #6b8e6b; color: white; text-decoration: none; border-radius: 5px; font-size: 0.95rem;'>RDKit</a>", unsafe_allow_html=True)
                with col2:
                    st.markdown(f"<a href='{info['ref2']}' target='_blank' style='display: block; text-align: center; padding: 0.5rem; background: #6b8e6b; color: white; text-decoration: none; border-radius: 5px; font-size: 0.95rem;'>Optuna</a>", unsafe_allow_html=True)
            elif 'ref' in info:
                # Cas normal: un seul lien
                st.markdown(f"<a href='{info['ref']}' target='_blank' style='display: block; text-align: center; padding: 0.5rem; background: #6b8e6b; color: white; text-decoration: none; border-radius: 5px; font-size: 0.95rem;'>En savoir plus</a>", unsafe_allow_html=True)
            st.markdown("</div>", unsafe_allow_html=True)



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
        'Rattus norvegicus': "Le rat norvégien (Rattus norvegicus), également appelé rat brun ou surmulot, est une espèce de rongeurs très utilisée en recherche biomédicale.\n\nC'est un organisme modèle important pour l'étude de la physiologie, de la pharmacologie et du comportement.\n\nSa proximité génétique avec l'homme en fait un modèle précieux pour la recherche sur les maladies humaines.",

        'Saccharomyces cerevisiae': "Saccharomyces cerevisiae, communément appelée levure de boulanger ou levure de bière, est un organisme eucaryote unicellulaire.\n\nC'est l'un des organismes modèles les plus étudiés en biologie cellulaire et moléculaire.\n\nSon génome a été le premier d'un eucaryote à être entièrement séquencé en 1996.",

        'Bos taurus': "Bos taurus est le nom scientifique du bœuf domestique ou bovin.\n\nBien que principalement élevé pour l'agriculture, il est également utilisé en recherche biomédicale, notamment pour l'étude du métabolisme, de la reproduction et de certaines maladies zoonotiques.",

        'Cavia porcellus': "Le cobaye (Cavia porcellus), également appelé cochon d'Inde, est un rongeur originaire d'Amérique du Sud.\n\nC'est un organisme modèle important en immunologie, en nutrition et en recherche sur les maladies infectieuses.\n\nContrairement à la plupart des mammifères, il ne peut pas synthétiser la vitamine C.",

        'Oryctolagus cuniculus': "Le lapin de garenne (Oryctolagus cuniculus) est utilisé comme organisme modèle en recherche biomédicale, notamment en immunologie, en toxicologie et en recherche cardiovasculaire.\n\nIl est particulièrement utile pour la production d'anticorps et l'étude de certaines maladies oculaires.",

        'Cricetulus griseus': "Le hamster chinois (Cricetulus griseus) est un petit rongeur originaire de Chine et de Mongolie.\n\nSes cellules ovariennes (cellules CHO) sont largement utilisées en biotechnologie pour la production de protéines thérapeutiques et en recherche sur la génétique cellulaire.",

        'Chlorocebus sabaeus': "Le singe vert (Chlorocebus sabaeus) est un primate originaire d'Afrique de l'Ouest.\n\nIl est utilisé comme organisme modèle en recherche biomédicale, notamment pour l'étude des maladies infectieuses, du VIH/SIDA et des troubles neurologiques en raison de sa proximité phylogénétique avec l'homme.",

        'Mus musculus': "La souris de laboratoire (Mus musculus) est l'organisme modèle mammifère le plus utilisé en recherche biomédicale.\n\nSon génome, très similaire à celui de l'homme, permet d'étudier de nombreuses maladies humaines.\n\nElle est utilisée dans pratiquement tous les domaines de la recherche biologique.",

        'Escherichia coli': "Escherichia coli (E. coli) est une bactérie intestinale commune qui constitue l'un des organismes modèles les plus importants en biologie moléculaire et en génétique microbienne.\n\nElle est largement utilisée en biotechnologie pour la production de protéines recombinantes et comme hôte pour le clonage génétique.",

        'Trypanosoma cruzi': "Trypanosoma cruzi est un protozoaire parasite responsable de la maladie de Chagas, une maladie tropicale négligée affectant des millions de personnes en Amérique latine.\n\nC'est un organisme modèle important pour l'étude des maladies parasitaires et le développement de nouveaux traitements antiparasitaires.",

        'Human immunodeficiency virus 1': "Le virus de l'immunodéficience humaine de type 1 (VIH-1) est le rétrovirus responsable du SIDA.\n\nIl est intensivement étudié en virologie et en immunologie pour comprendre les mécanismes d'infection virale et développer des thérapies antivirales et des vaccins.",

        'Mycobacterium tuberculosis': "Mycobacterium tuberculosis est la bactérie responsable de la tuberculose, l'une des maladies infectieuses les plus mortelles au monde.\n\nC'est un organisme modèle crucial pour la recherche sur les maladies infectieuses, le développement d'antibiotiques et l'étude des mécanismes de résistance aux médicaments.",

        'Drosophila': "La drosophile (Drosophila melanogaster), ou mouche du vinaigre, est un organisme modèle fondamental en génétique et en biologie du développement.\n\nSon cycle de vie court, sa facilité de culture et son génome bien caractérisé en font un outil précieux pour étudier l'hérédité, le développement embryonnaire et les maladies neurodégénératives.",

        'Schistosoma mansoni': "Schistosoma mansoni est un ver parasite responsable de la schistosomiase, une maladie tropicale affectant des millions de personnes.\n\nC'est un organisme modèle important pour l'étude des helminthiases et le développement de médicaments antiparasitaires.",

        'Caenorhabditis elegans': "Caenorhabditis elegans est un nématode transparent d'environ 1 mm de long.\n\nC'est un organisme modèle majeur en biologie du développement, en génétique et en neurobiologie.\n\nSon système nerveux simple mais complet, composé de seulement 302 neurones, est entièrement cartographié.",

        'Bacillus subtilis': "Bacillus subtilis est une bactérie Gram-positive couramment trouvée dans le sol.\n\nC'est un organisme modèle important en microbiologie pour l'étude de la formation des spores, de la différenciation cellulaire et de la production d'enzymes industrielles.",

        'Leishmania mexicana': "Leishmania mexicana est un protozoaire parasite responsable de la leishmaniose cutanée.\n\nC'est un organisme modèle pour l'étude des maladies parasitaires tropicales et le développement de nouveaux traitements contre les leishmanioses.",

        'Trypanosoma brucei TREU927': "Trypanosoma brucei est un protozoaire parasite responsable de la maladie du sommeil en Afrique.\n\nLa souche TREU927 est utilisée comme organisme modèle de référence pour l'étude de la biologie des trypanosomes et le développement de médicaments antiparasitaires."
    }

    # Retourner l'information ou un message par défaut
    return organism_info.get(organism_name, f"Organisme modèle utilisé en recherche biomédicale pour l'étude de diverses propriétés biologiques et le développement de nouveaux composés thérapeutiques.")

def call_api(smiles: str, organism: str = "Homo sapiens") -> Dict[str, Any]:
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

        # Paramètres de la requête
        payload = {
            "smiles": smiles,
            "organism": organism
        }

        # Appel API
        response = requests.post(url, json=payload, timeout=30)

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
    st.markdown("<p style='text-align: center; color: #b0b0b0; font-size: 1rem;'>Propriétés disponibles : Stress Oxydatif, Métabolisme Énergétique, Mort Cellulaire et Signalisation Cellulaire</p>", unsafe_allow_html=True)
    st.markdown("---")

    # Sidebar - Informations et paramètres
    with st.sidebar:
        # Exemples de SMILES
        st.markdown("### 📋 Exemples de SMILES")
        examples = {
            "Aspirine": "CC(=O)Oc1ccccc1C(=O)O",
            "Caféine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "Pénicilline G": "CC1(C(N2C(S1)C(C2=O)NC(=O)Cc3ccccc3)C(=O)O)C"
        }

        for name, smiles in examples.items():
            if st.button(f"{name}", use_container_width=True):
                st.session_state['example_smiles'] = smiles

        st.markdown("---")
        st.markdown("### 📖 À propos")
        st.markdown("""
        **BioGNN** est un modèle prédictif d'activité biologique de molécules sur des organismes modèles.\n
        Projet réalisé dans le cadre des projets du Bootcamp Data Science & IA de Le Wagon.
        """)
        st.markdown("---")
        # Bouton pour afficher l'équipe
        if st.button("👥 About the team", use_container_width=True):
            show_team_modal()
        if st.button("🧠 About the model", use_container_width=True):
            show_model_info()

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

    # Gérer les résultats dans session_state pour éviter les rechargements
    if predict_button and smiles_input:
        is_valid, _ = validate_smiles(smiles_input)

        if not is_valid:
            st.error("⚠️ Veuillez entrer un SMILES valide avant de prédire")
        else:
            # Appel à l'API

            with st.spinner("🔄 Analyse en cours..."):
                result = call_api(smiles_input, selected_organism)

            # Stocker les résultats dans session_state
            if result["success"]:
                st.session_state['last_prediction'] = {
                    'result': result,
                    'smiles': smiles_input,
                    'organism': selected_organism
                }
            else:
                st.session_state['last_prediction'] = None
                st.error(f"❌ {result['error']}")

    # Afficher les résultats s'ils existent dans session_state
    if 'last_prediction' in st.session_state and st.session_state['last_prediction'] is not None:
        pred = st.session_state['last_prediction']
        result = pred['result']
        smiles_input = pred['smiles']
        selected_organism = pred['organism']

        # Affichage de la prédiction

        if result["success"]:
            data = result["data"]

            # ===============================
            # 🔬 Vérification toxicité si organisme supporté
            # ===============================
            selected_organism_tox = selected_organism

            if selected_organism_tox in ORGANISMS_TOXICITY_MAPPING.values():
                tox_payload = {
                    "smiles": smiles_input,
                    "organism": selected_organism_tox
                }

                try:
                    tox_response = requests.post(
                        f"{BASE_URI}predict_tox",
                        json=tox_payload,
                        timeout=15
                    )

                    if tox_response.status_code == 200:
                        tox_data = tox_response.json()

                        prob_toxicity = tox_data.get("prob_toxicity", 0.0)
                        toxic = tox_data.get("toxic", False)

                        if toxic:
                            st.error(
                                f"La molécule d'intérêt a une probabilité de {prob_toxicity:.3f} "
                                f"d'être toxique pour l'organisme étudié. Soyez attentifs aux dosages."
                            )

                    else:
                        st.info("ℹ️ Analyse de toxicité non disponible pour cet organisme")

                except Exception:
                    st.info("ℹ️ Impossible de récupérer la prédiction de toxicité")

            properties = data.get("properties", {})

            # Générer dynamiquement le summary basé sur un threshold
            threshold = 0.8
            promising_properties = [prop for prop, score in properties.items() if score >= threshold]

            if promising_properties:
                summary = f"Candidat prometteur pour : {', '.join(promising_properties)}"
            else:
                summary = f"Aucune propriété fortement prédite (probabilités < {threshold:.2f})"

            st.markdown(
                f'<p class="prediction-text" style="font-size:1.5rem;">{summary}</p>',
                unsafe_allow_html=True
            )

            st.markdown("#### 🧪 Scores par propriété biologique")

            # Affichage en grille 2 colonnes
            props_list = list(properties.items())
            cols_per_row = 2
            for i in range(0, len(props_list), cols_per_row):
                row_items = props_list[i:i+cols_per_row]
                cols = st.columns(len(row_items))
                for col, (prop, score) in zip(cols, row_items):
                    if score >= 0.7:
                        color = "#b8e986"
                    elif score >= 0.4:
                        color = "#f0d264"
                    else:
                        color = "#c0c0c0"

                    with col:
                        st.markdown(
                            f"""
                            <div class="result-card">
                                <strong>{prop}</strong>
                                <div style="margin-top:0.4rem;">
                                    Probabilité d'activité :
                                    <strong style="color:{color};">
                                        {score:.3f}
                                    </strong>
                                </div>
                            </div>
                            """,
                            unsafe_allow_html=True
                        )

            st.markdown('</div>', unsafe_allow_html=True)

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

            # Affichage : 2 colonnes verticales alignées - Organisme + Propriétés

            def get_property_description(property_name: str) -> str:
                """
                Retourne la description détaillée d'une propriété biologique prédite par le modèle BioGNN.
                """
                descriptions = {
                "Stress Oxydatif": "Prédiction de l'activité antioxydante et de la modulation du stress oxydatif. Évalue la capacité de la molécule à agir comme antioxydant direct ou indirect.\n\nUn score élevé suggère un potentiel pour des applications nutraceutiques, cosmétiques ou agrochimiques.\n\nNote : la corrélation entre activité in vitro et efficacité in vivo reste un défi majeur dans ce domaine.",

                "Mort Cellulaire": "Prédiction de l'induction de mort cellulaire par différents mécanismes : apoptose, autophagie ou ferroptose. Propriété pertinente pour la découverte de composés anticancéreux, où l'induction sélective dans les cellules tumorales est recherchée.\n\nUn score élevé indique un potentiel cytotoxique significatif.\n\nNote : le mécanisme précis de mort cellulaire et le contexte cellulaire influencent fortement l'activité biologique réelle.",

                "Métabolisme Énergétique": "Prédiction de l'activité sur le métabolisme énergétique cellulaire, incluant la glycolyse et la fonction mitochondriale. Cette propriété capture le potentiel de modulation du métabolisme glucidique et lipidique, avec des applications dans le diabète, l'obésité et le métabolisme tumoral (effet Warburg).\n\nUn score élevé suggère une interaction avec les voies métaboliques centrales.\n\nAttention : la toxicité mitochondriale constitue un effet indésirable majeur à considérer, et l'effet peut varier significativement selon le contexte tissulaire (foie, muscle, tissu adipeux).",

                "Signalisation Cellulaire": "Prédiction de l'activité sur les grandes voies de signalisation intracellulaires, en évaluant le potentiel d'inhibition de kinases ou de modulation de GPCR.\n\nUn score élevé indique une probable interaction avec ces cascades de signalisation qui régulent la prolifération, la différenciation et la réponse immunitaire."
                }
                return descriptions.get(property_name, "Description non disponible pour cette propriété.")

            # 2 colonnes de même taille
            col_organism, col_properties = st.columns([1, 1])

            # Colonne 1: Organisme
            with col_organism:
                st.markdown(f"""
                <div class="prediction-card" style="text-align: center; min-height: 440px;">
                <p class="prediction-label">ORGANISME SÉLECTIONNÉ</p>
                <p class="prediction-text" style="font-size: 1.8rem;">{selected_organism}</p>
                <div class="info-box" style="margin-top: 1rem; text-align: left;">
                    <p style="margin: 0.5rem 0; color: #e0e0e0;">
                        {get_organism_info(selected_organism)}
                    </p>
                </div>
                </div>
                """, unsafe_allow_html=True)

            # Colonne 2: Propriétés avec sélecteur
            with col_properties:
                # Utiliser un container Streamlit avec un style CSS global
                st.markdown("""
                <style>
                /* Cibler spécifiquement cette colonne */
                div[data-testid="column"]:has(#properties-title) {
                    background: linear-gradient(135deg, #4a5d4e 0%, #5a6d5e 100%) !important;
                    padding: 1.5rem !important;
                    border-radius: 15px !important;
                    box-shadow: 0 4px 6px rgba(0, 0, 0, 0.3) !important;
                    border: 2px solid #6b7d6b !important;
                    min-height: 400px !important;
                }
                </style>
                </div>
                """, unsafe_allow_html=True)

                # Sélecteur de propriété
                prop_names = list(properties.keys())
                prop_options = [f"{name} (Score: {properties[name]:.2f})" for name in prop_names]

                # Initialiser l'index de la propriété sélectionnée dans session_state si nécessaire
                if 'selected_prop_index' not in st.session_state:
                    st.session_state.selected_prop_index = 0

                selected_prop_display = st.selectbox(
                    "Sélectionnez une propriété:",
                    prop_options,
                    index=st.session_state.selected_prop_index,
                    key=f"property_selector_{smiles_input}_{selected_organism}",
                    label_visibility="visible"
                )

                # Mettre à jour l'index dans session_state
                st.session_state.selected_prop_index = prop_options.index(selected_prop_display)

                # Extraire le nom de la propriété sélectionnée
                selected_prop_name = prop_names[prop_options.index(selected_prop_display)]
                selected_prop_score = properties[selected_prop_name]

                # Déterminer la couleur selon le score
                if selected_prop_score >= 0.7:
                    score_color = "#b8e986"
                    score_label = "Fort"
                elif selected_prop_score >= 0.4:
                    score_color = "#f0d264"
                    score_label = "Modéré"
                else:
                    score_color = "#c0c0c0"
                    score_label = "Faible"

                # Afficher la description de la propriété sélectionnée
                st.markdown(f"""
                <div class="info-box" style="margin-top: 0.1rem;">
                    <h4 style="color: #b8e986; margin-top: 0; margin-bottom: 0.3rem;">{selected_prop_name}</h4>
                    <p style="margin: 0.8rem 0; color: #e0e0e0; line-height: 1.6;">
                        {get_property_description(selected_prop_name)}
                    </p>
                    <div style="margin-top: 1rem; padding: 0.8rem; background-color: #3d4a3e; border-radius: 8px;">
                        <p style="margin: 0; color: #d4d4d4;">
                            <strong>Score prédit:</strong>
                            <span style="color: {score_color}; font-size: 1.2rem; font-weight: bold;">
                                {selected_prop_score:.2f}
                            </span>
                            <span style="color: {score_color}; margin-left: 0.5rem;">
                                ({score_label})
                            </span>
                        </p>
                    </div>
                </div>
                """, unsafe_allow_html=True)

        else:
            st.error(f"❌ {result['error']}")

    # Message si le bouton est cliqué sans SMILES valide
    elif predict_button:
        st.warning("⚠️ Veuillez entrer un SMILES valide")

    # Footer
    st.markdown("""
    <div style='text-align: center; color: #808080; padding: 2rem;'>
        <p>🧬 BioGNN - Deep Learning for Molecular Discovery 🧬</p>
        <p>Les résultats ne sont que prédictifs et ne visent qu'à aiguiller les décisions en amont d'essais expérimentaux.</p>
    </div>
    """, unsafe_allow_html=True)

# ============================================================================
# POINT D'ENTRÉE
# ============================================================================

if __name__ == "__main__":
    main()
