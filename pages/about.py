import streamlit as st

# Page principale
def main():
    st.title("Ma superbe application Streamlit")
    # Votre contenu ici

# Page de retour
def retour_page_principale():
    st.title("Retour à la page principale")
    st.write("Vous êtes revenu à la page principale.")

# Vérifier si le bouton "Retour" est cliqué
if st.button("Retour"):
    # Revenir à la page principale
    st.experimental_rerun()

# Vérifier quelle page afficher
if __name__ == "__main__":
    if st.session_state.page == "retour":
        retour_page_principale()
    else:
        main()