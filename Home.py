import streamlit as st
from openai import OpenAI
from openai.types.beta.assistant_stream_event import ThreadMessageDelta
from openai.types.beta.threads.text_delta_block import TextDeltaBlock

# Load API keys and assistant ID from Streamlit secrets
OPENAI_API_KEY = st.secrets["OPENAI_API_KEY"]
ASSISTANT_ID = st.secrets["ASSISTANT_ID"]

# Initialise the OpenAI client and retrieve the assistant
client = OpenAI(api_key=OPENAI_API_KEY)
assistant = client.beta.assistants.retrieve(assistant_id=ASSISTANT_ID)

# Initialise session state to store conversation history and track run state
if "chat_history" not in st.session_state:
    st.session_state.chat_history = []
if "run_active" not in st.session_state:
    st.session_state.run_active = False

# Page title and introduction
st.title("Andreeff's Lab Publications Assistant")
st.write("Welcome! Ask me any question about Dr. Andreeff's Lab publications and I'll try to help.")

# Display chat history
for message in st.session_state.chat_history:
    with st.chat_message(message["role"]):
        st.markdown(message["content"])

# Text input for user queries
user_query = st.chat_input("Ask me a question about the publications:", key="query_input")

if user_query:
    if st.session_state.run_active:
        # Notify user if the assistant is still processing a previous query
        st.warning("The assistant is currently processing your previous question. Please wait.")
    else:
        # Indicate that a processing run has started
        st.session_state.run_active = True

        # Create a new chat thread if it doesn't already exist
        if "thread_id" not in st.session_state:
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
                assistant_id=ASSISTANT_ID,
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
