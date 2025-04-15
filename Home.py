import streamlit as st
from openai import OpenAI
from openai.types.beta.assistant_stream_event import ThreadMessageDelta
from openai.types.beta.threads.text_delta_block import TextDeltaBlock
st.markdown("""
<style>
div[data-testid="stToolbar"] {display: none;}
</style>
""", unsafe_allow_html=True)

# Load API keys and assistant ID from Streamlit secrets
OPENAI_API_KEY = st.secrets["OPENAI_API_KEY"]
client = OpenAI(api_key=OPENAI_API_KEY)
ASSISTANT_ID_ASH = "asst_O5eLrcWgvxDTSWz1kJks3blQ"  # ASH 2024 abstracts
ASSISTANT_ID_MHT = "asst_KTkXrKD646vjDo2W9IrzfnK3"  # MHT publications
ASSISTANT_ID_HACKATHON = "asst_pFOiqiKqq0xombtZJLKQpXJ4"  # Hackathon

# Initialize session state to store conversation history and track run state
if "chat_history" not in st.session_state:
    st.session_state.chat_history = []
if "run_active" not in st.session_state:
    st.session_state.run_active = False
if "assistant_id" not in st.session_state:
    st.session_state.assistant_id = ASSISTANT_ID_HACKATHON  # Default to MHT publications
if "thread_id" not in st.session_state:
    st.session_state.thread_id = None

# Function to reset the chat history and thread
def reset_chat():
    st.session_state.chat_history = []
    st.session_state.thread_id = None
    st.session_state.run_active = False

# Sidebar for user instructions and assistant selection
with st.sidebar:
    st.markdown(
        """
        ### How to Use
        - **MDACC 2025 Hackathon**: Ask questions about Hackathon papers and jupyter notebook.
        - The assistant will provide responses based on the selected data source.

        ### Disclaimer
        - The AI might make mistakes; always verify the information provided.
        """,
        unsafe_allow_html=True,
    )
    # - Choose an assistant to interact with from the dropdown menu.

    # - **ASH 2024 Abstracts**: Ask questions about ASH 2024 conference abstracts.
    # - **MHT Publications**: Ask questions about Dr. Andreeff's Lab publications.

    # Dropdown for selecting the assistant
    assistant_choice = st.selectbox(
        "Select an Assistant",
        ("ASH 2024 Abstracts", "MHT Publications"), # "MDACC 2025 Hackathon"
        index=0
    )

    # Set assistant ID based on user choice
    if assistant_choice == "MDACC 2025 Hackathon":
        st.session_state.assistant_id = ASSISTANT_ID_HACKATHON
    elif assistant_choice == "ASH 2024 Abstracts":
        st.session_state.assistant_id = ASSISTANT_ID_ASH
    else:
        st.session_state.assistant_id = ASSISTANT_ID_MHT

    # Add a reset button
    if st.button("Reset Chat"):
        reset_chat()
        st.rerun() # Rerun the app to clear the chat display

# Page title and introduction
st.title(f"{assistant_choice} Assistant")

# Display chat history
for message in st.session_state.chat_history:
    with st.chat_message(message["role"]):
        st.markdown(message["content"])

# Text input for user queries
user_query = st.chat_input("Ask me a question:", key="query_input")

if user_query:
    if st.session_state.run_active:
        # Notify user if the assistant is still processing a previous query
        st.warning("The assistant is currently processing your previous question. Please wait.")
    else:
        # Indicate that a processing run has started
        st.session_state.run_active = True

        # Create a new chat thread if it doesn't already exist OR if the chat was reset
        if "thread_id" not in st.session_state or st.session_state.thread_id is None:
            thread = client.beta.threads.create()
            st.session_state.thread_id = thread.id

        # Display user query in chat history
        with st.chat_message("user"):
            st.markdown(user_query)

        # Log user query to the history
        st.session_state.chat_history.append({"role": "user", "content": user_query})

        # Submit user query to the thread
        client.beta.threads.messages.create(
            thread_id=st.session_state.thread_id,
            role="user",
            content=user_query
        )

        # Handle streaming of assistant responses
        with st.chat_message("assistant"):
            stream = client.beta.threads.runs.create(
                thread_id=st.session_state.thread_id,
                assistant_id=st.session_state.assistant_id,  # Use selected assistant ID
                stream=True
            )

            # Initialize an empty container for the assistant's responses
            assistant_reply_box = st.empty()
            assistant_reply = ""

            # Iterate through the streaming events
            for event in stream:
                if isinstance(event, ThreadMessageDelta):
                    assistant_reply += event.data.delta.content[0].text.value
                    assistant_reply_box.markdown(assistant_reply)

            # Update chat history and reset run state after completion
            st.session_state.chat_history.append({"role": "assistant", "content": assistant_reply})
            st.session_state.run_active = False
