import streamlit as st

# Titre de l'application
st.title("La page about pour tester")

# Sous-titre
st.markdown("## Interface utilisateur interactive")

# Zone de texte
user_input = st.text_input("test de la deuxième page")

# Bouton
button_clicked = st.button("Cliquez ici")
