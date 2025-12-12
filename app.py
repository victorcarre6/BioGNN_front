"""
BioGNN Immunity - Application Streamlit
Interface utilisateur pour la prédiction de propriétés biologiques de molécules
"""

import os
import streamlit as st
import requests
from typing import Tuple, Dict, Optional, Any

# Import conditionnel de RDKit (peut ne pas être disponible sur certaines plateformes)
try:
    from rdkit import Chem
    from rdkit.Chem import Draw, Descriptors
    from PIL import Image
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    st.warning("⚠️ RDKit n'est pas disponible. La visualisation moléculaire sera limitée.")

# ============================================================================
# CONFIGURATION
# ============================================================================

st.set_page_config(
    page_title="BioGNN Immunity",
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

def smiles_to_image(smiles: str, size=(300, 300)) -> Optional[Any]:
    """
    Convertit un SMILES en image de molécule
    """
    if not RDKIT_AVAILABLE:
        return None

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        img = Draw.MolToImage(mol, size=size)
        return img
    except Exception as e:
        st.error(f"Erreur lors de la génération de l'image: {e}")
        return None

def get_molecule_properties(smiles: str) -> Dict[str, Any]:
    """
    Calcule les propriétés de base d'une molécule
    """
    if not RDKIT_AVAILABLE:
        return {
            "SMILES": smiles,
            "Note": "RDKit non disponible - propriétés limitées"
        }

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {}

        return {
            "Masse moléculaire": f"{Descriptors.MolWt(mol):.2f} g/mol",
            "LogP": f"{Descriptors.MolLogP(mol):.2f}",
            "Nombre d'atomes": mol.GetNumAtoms(),
            "Nombre de liaisons": mol.GetNumBonds(),
            "Donneurs H": Descriptors.NumHDonors(mol),
            "Accepteurs H": Descriptors.NumHAcceptors(mol),
            "Cycles aromatiques": Descriptors.NumAromaticRings(mol),
        }
    except Exception as e:
        st.error(f"Erreur calcul propriétés: {e}")
        return {}

def call_api(smiles: str, property_name: str) -> Dict[str, Any]:
    """
    Appelle l'API GCP pour obtenir une prédiction

    TODO: Adapter les paramètres selon votre API
    """
    try:
        url = BASE_URI + "predict"

        # Paramètres de la requête - ADAPTER SELON VOTRE API
        params = {
            "smiles": smiles,
            "property": property_name
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
    st.markdown("<h1 style='text-align: center;'>🧬 BioGNN Immunity</h1>", unsafe_allow_html=True)
    st.markdown("<p style='text-align: center; color: #b0b0b0; font-size: 1.2rem;'>Prédiction de propriétés biologiques par Graph Neural Networks</p>", unsafe_allow_html=True)

    st.markdown("---")

    # Sidebar - Informations et paramètres
    with st.sidebar:
        st.markdown("### ⚙️ Configuration")

        # Afficher l'URL de l'API (pour debug)
        with st.expander("🔗 API Configuration"):
            st.code(BASE_URI, language="text")

            # Test de connexion
            if st.button("🔍 Tester la connexion API"):
                try:
                    health_url = BASE_URI.rstrip('/') + '/health'
                    response = requests.get(health_url, timeout=5)
                    if response.status_code == 200:
                        st.success("✅ API accessible")
                    else:
                        st.error(f"⚠️ Status: {response.status_code}")
                except Exception as e:
                    st.error(f"❌ Erreur: {str(e)}")

        st.markdown("---")

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
        **BioGNN Immunity** utilise des Graph Neural Networks pour prédire
        les propriétés biologiques de molécules à partir de leur structure SMILES.
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
        st.markdown("#### 🎯 Propriété Cible")

        # Sélection de propriété
        properties = [
            "Activité antimicrobienne",
            "Stress oxydatif & défenses",
            "Cycle cellulaire & prolifération",
            "Mort cellulaire",
            "Inflammation & immunité",
            "Signalisation cellulaire",
            "Intégrité génomique",
            "Métabolisme énergétique",
            "Homéostasie tissulaire",
            "Fonctions spécifiques d'organes"
        ]

        selected_property = st.selectbox(
            "Sélectionnez la propriété à prédire:",
            properties,
            label_visibility="collapsed"
        )

    # Bouton de prédiction
    st.markdown("<br>", unsafe_allow_html=True)
    col_btn1, col_btn2, col_btn3 = st.columns([1, 2, 1])
    with col_btn2:
        predict_button = st.button("🔮 PRÉDIRE", use_container_width=True)

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
            st.markdown("### 🔬 Affichage de la molécule")

            col_mol1, col_mol2 = st.columns([1, 1])

            with col_mol1:
                st.markdown('<div class="molecule-container">', unsafe_allow_html=True)
                if RDKIT_AVAILABLE:
                    mol_img = smiles_to_image(smiles_input, size=(400, 400))
                    if mol_img:
                        st.image(mol_img, use_container_width=True)
                    else:
                        st.info("🧪 Impossible de générer l'image de la molécule")
                else:
                    st.info("🧪 **Visualisation moléculaire non disponible**\n\nRDKit n'est pas installé. La molécule sera traitée par l'API backend.")
                    st.code(smiles_input, language="text")
                st.markdown('</div>', unsafe_allow_html=True)

            with col_mol2:
                st.markdown("#### 📊 Propriétés moléculaires")
                props = get_molecule_properties(smiles_input)

                if props:
                    for key, value in props.items():
                        st.markdown(f"""
                        <div class="info-box">
                            <strong>{key}:</strong> {value}
                        </div>
                        """, unsafe_allow_html=True)

            st.markdown("---")

            # Appel à l'API
            st.markdown("### 🎯 Prédiction")

            with st.spinner("🔄 Analyse en cours..."):
                result = call_api(smiles_input, selected_property)

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
        <p>🧬 BioGNN Immunity - Propulsé par Graph Neural Networks</p>
        <p style='font-size: 0.9rem;'>Développé avec ❤️ pour l'immunité computationnelle</p>
    </div>
    """, unsafe_allow_html=True)

# ============================================================================
# POINT D'ENTRÉE
# ============================================================================

if __name__ == "__main__":
    main()
