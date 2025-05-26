
try:
    # Read as binary
    with open("E:\\Development\\PyMISES\\docs\\status.md", 'rb') as f:
        binary_content = f.read()
    
    # Try to decode with replacement
    text_content = binary_content.decode('utf-8', errors='replace')
    
    # Write decoded content
    with open("E:\\Development\\PyMISES\\docs\\status_clean.txt", 'w', encoding='utf-8') as f:
        f.write(text_content)
    
    print("Successfully created status_clean.txt")
    
except Exception as e:
    print(f"Error: {e}")
