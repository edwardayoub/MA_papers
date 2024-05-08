import streamlit as st
from openai import OpenAI

OPENAI_API_KEY = st.secrets["OPENAI_API_KEY"]
ASSISTANT_ID = st.secrets["ASSISTANT_ID"]

# Initialise the OpenAI client, and retrieve the assistant
client = OpenAI(api_key=OPENAI_API_KEY)
assistant = client.beta.assistants.retrieve(assistant_id=ASSISTANT_ID)

# Initialise session state to store conversation history and run state
if "chat_history" not in st.session_state:
    st.session_state.chat_history = []
if "run_active" not in st.session_state:
    st.session_state.run_active = False

# Title
st.title("MA publications chatbot")

# Display messages in chat history
for message in st.session_state.chat_history:
    with st.chat_message(message["role"]):
        st.markdown(message["content"])

# Textbox and streaming process
user_query = st.chat_input("Ask me a question", key="query_input")

if user_query:
    if st.session_state.run_active:
        # Display a message to the user indicating that the assistant is still processing
        st.warning("The assistant is still processing your previous question. Please wait.")
    else:
        # Update state to indicate a run is active
        st.session_state.run_active = True

        # Create a new thread if it does not exist
        if "thread_id" not in st.session_state:
            thread = client.beta.threads.create()
            st.session_state.thread_id = thread.id

        # Display the user's query
        with st.chat_message("user"):
            st.markdown(user_query)

        # Store the user's query into the history
        st.session_state.chat_history.append({"role": "user", "content": user_query})

        # Add user query to the thread
        client.beta.threads.messages.create(
            thread_id=st.session_state.thread_id,
            role="user",
            content=user_query
        )

        # Stream the assistant's reply
        with st.chat_message("assistant"):
            stream = client.beta.threads.runs.create(
                thread_id=st.session_state.thread_id,
                assistant_id=ASSISTANT_ID,
                stream=True
            )

            # Empty container to display the assistant's reply
            assistant_reply_box = st.empty()

            # A blank string to store the assistant's reply
            assistant_reply = ""

            # Iterate through the stream 
            for event in stream:
                if isinstance(event, ThreadMessageDelta):
                    # Process text delta blocks
                    assistant_reply += event.data.delta.content[0].text.value
                    assistant_reply_box.markdown(assistant_reply)

            # Once the stream is over, update chat history and reset run_active state
            st.session_state.chat_history.append({"role": "assistant", "content": assistant_reply})
            st.session_state.run_active = False
