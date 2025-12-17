# BioGNN

---

Deep learning model to predict biological activity of organic molecules. The project leverages Graph Neural Networks (GNN) to exploit the natural representation of molecules as graphs (atoms = nodes, bonds = edges), with a pipeline that include data collection from ChemBL and feature extraction from RDKit.

This project is part of the final project for the Le Wagon Data Science & AI bootcamp.

---

## 🚀 Configuration et Déploiement

### Configuration de l'API

L'application peut fonctionner avec différentes configurations d'API selon l'environnement:

#### 1. **Production (Streamlit Cloud)**
Par défaut, l'application utilise l'URL GCP définie dans le code:
```
https://biognn-api-271034908172.europe-west1.run.app/
```

Pour Streamlit Cloud, copiez dans Settings > Secrets:
```toml
cloud_api_uri = "https://biognn-api-271034908172.europe-west1.run.app/"
```

#### 2. **Développement Local avec API Cloud**
Aucune configuration nécessaire. L'app utilisera automatiquement l'API cloud par défaut.

#### 3. **Développement Local avec API Locale**

**Option A: Via secrets.toml**
1. Copiez `.streamlit/secrets.toml.example` vers `.streamlit/secrets.toml`
2. Décommentez la ligne:
   ```toml
   use_local_api = true
   ```
3. Vérifiez que `local_api_uri` pointe vers votre API locale:
   ```toml
   local_api_uri = "http://localhost:8000"
   ```

**Option B: Via variable d'environnement (priorité maximale)**
```bash
export API_URL="http://localhost:8000/"
streamlit run app_update.py
```

#### 4. **Docker Local**
1. Décommentez dans `.streamlit/secrets.toml`:
   ```toml
   use_docker_api = true
   ```
2. Vérifiez `local_docker_uri`:
   ```toml
   local_docker_uri = "http://localhost:8080"
   ```

### Ordre de Priorité des URLs

L'application sélectionne l'URL de l'API dans l'ordre suivant:
1. Variable d'environnement `API_URL`
2. `local_api_uri` (si `use_local_api = true`)
3. `local_docker_uri` (si `use_docker_api = true`)
4. `cloud_api_uri` (si défini dans secrets)
5. URL GCP par défaut (hardcodée)

### Vérifier la Configuration

L'URL utilisée est affichée dans la sidebar sous "🔧 Configuration API". Ouvrez cet expander pour voir:
- L'URL actuelle de l'API
- L'environnement détecté (Local ou Production)

---

