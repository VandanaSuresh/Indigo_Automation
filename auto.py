import os
import time
import re
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

# Updated Paths
input_folder_path = r"C:\Users\91934\Desktop\TIGS\Automation_new\INPUT"
wild_type_file_path = r"C:\Users\91934\Desktop\TIGS\Automation_new\WILD_TYPE\wildtype.ab1"
download_dir = r"C:\Users\91934\Desktop\TIGS\Automation_new\NEW"

# Ensure input folder exists
if not os.path.exists(input_folder_path):
    print(f"Error: Input folder '{input_folder_path}' does not exist.")
    exit(1)

# Set up ChromeDriver
service = Service(r"C:\New folder\chromedriver-win64\chromedriver.exe")
options = webdriver.ChromeOptions()
prefs = {"download.default_directory": download_dir}
options.add_experimental_option("prefs", prefs)

driver = webdriver.Chrome(service=service, options=options)

# Function to highlight multiple PAM sequences in HTML
def highlight_pam_sequence(html_content, pam_sequences=None):
    if pam_sequences is None:
        pam_sequences = ["CGGAGGAC"]  # Default PAM sequence

    for pam in pam_sequences:
        pattern = re.compile(re.escape(pam), re.IGNORECASE)
        html_content = pattern.sub(
            lambda m: f'<span style="background-color: yellow; font-weight: bold;">{m.group(0)}</span>',
            html_content
        )
    return html_content

# Function to process a single file
def process_input_file(input_file_path, pam_sequences):
    input_file_base_name = os.path.splitext(os.path.basename(input_file_path))[0]
    driver.get("https://www.gear-genomics.com/indigo/")
    wait = WebDriverWait(driver, 20)
    
    try:
        input_file_upload = wait.until(EC.presence_of_element_located((By.ID, "inputFile")))
        input_file_upload.send_keys(input_file_path)
        print(f"Uploaded {input_file_base_name}")
        
        wild_type_link = wait.until(EC.element_to_be_clickable((By.ID, "target-chromatogram-tab")))
        driver.execute_script("arguments[0].scrollIntoView();", wild_type_link)
        time.sleep(1)
        driver.execute_script("arguments[0].click();", wild_type_link)
        
        wild_type_file_input = wait.until(EC.presence_of_element_located((By.ID, "targetFileChromatogram")))
        wild_type_file_input.send_keys(wild_type_file_path)
        
        submit_button = wait.until(EC.element_to_be_clickable((By.ID, "btn-submit")))
        driver.execute_script("arguments[0].scrollIntoView();", submit_button)
        time.sleep(1)
        driver.execute_script("arguments[0].click();", submit_button)
        
        time.sleep(5)
        try:
            result_download_link = wait.until(EC.presence_of_element_located((By.LINK_TEXT, "Download HTML")))
            driver.execute_script("arguments[0].click();", result_download_link)
            time.sleep(10)
        except:
            results_html = driver.page_source
            highlighted_html = highlight_pam_sequence(results_html, pam_sequences=pam_sequences)
            result_html_file_path = os.path.join(download_dir, f"{input_file_base_name}_results_highlighted.html")
            with open(result_html_file_path, "w", encoding="utf-8") as html_file:
                html_file.write(highlighted_html)
            print(f"Saved results for {input_file_base_name}")
    except Exception as e:
        print(f"Error processing {input_file_base_name}: {e}")

# List of PAM sequences to highlight
pam_list = ["AGGAGGAC"]  # Add more sequences here as needed

# Process all .ab1 files
for file_name in os.listdir(input_folder_path):
    if file_name.endswith(".ab1"):
        input_file_path = os.path.join(input_folder_path, file_name)
        process_input_file(input_file_path, pam_sequences=pam_list)

driver.quit()
