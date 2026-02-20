# created by Vandana Suresh
# date - january 2025

import os
import time
import re
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from datetime import datetime
from selenium.common.exceptions import InvalidSessionIdException

# Updated Paths
input_folder_path = r"C:\Users\91934\Desktop\TIGS\Automation_new\INPUT\RG317"
wild_type_file_path = r"C:\Users\91934\Desktop\TIGS\Automation_new\WILD_TYPE\wildtype.ab1"
download_dir = r"C:\Users\91934\Desktop\TIGS\Automation_new\NEW"
log_file_path = os.path.join(download_dir, "analysis_log.txt")

# Ensure input folder exists
if not os.path.exists(input_folder_path):
    print(f"Error: Input folder '{input_folder_path}' does not exist.")
    exit(1)

# Function to initialize WebDriver
def init_driver():
    service = Service(r"C:\New_folder\chromedriver-win64\chromedriver-win64\chromedriver.exe")
    service.log_path = os.devnull  # Suppress ChromeDriver logs
    options = webdriver.ChromeOptions()
    prefs = {"download.default_directory": download_dir}
    options.add_experimental_option("prefs", prefs)
    options.add_argument("--log-level=3")  # Minimize browser logging
    options.add_experimental_option("excludeSwitches", ["enable-logging"])  # Suppress DevTools errors
    # options.add_argument("--headless=new")  # Uncomment for headless mode (test first)
    return webdriver.Chrome(service=service, options=options)

# Initialize WebDriver
driver = init_driver()

# Function to highlight multiple target gRNA sequences in HTML with different colors
def highlight_pam_sequence(html_content, grna_sequences=None):
    if grna_sequences is None:
        grna_sequences = ["CAGCAGCTGG","TCCAACCAGG"]  # Default gRNA sequence
    colors = ["cyan","yellow"]
    for i, grna_sequence in enumerate(grna_sequences):
        color = colors[i % len(colors)]
        pattern = re.compile(re.escape(grna_sequence), re.IGNORECASE)
        html_content = pattern.sub(
            lambda m: f'<span style="background-color: {color}; font-weight: bold;">{m.group(0)}</span>',
            html_content
        )
    return html_content

# Function to process a single file
def process_input_file(input_file_path, grna_sequences, driver):
    input_file_base_name = os.path.splitext(os.path.basename(input_file_path))[0]
    try:
        driver.get("https://www.gear-genomics.com/indigo/")
        wait = WebDriverWait(driver, 20)
       
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
            highlighted_html = highlight_pam_sequence(results_html, grna_sequences=grna_sequences)
            result_html_file_path = os.path.join(download_dir, f"{input_file_base_name}_results_highlighted.html")
            with open(result_html_file_path, "w", encoding="utf-8") as html_file:
                html_file.write(highlighted_html)
            print(f"Saved results for {input_file_base_name}")
        return True  # Success
    except Exception as e:
        print(f"Error processing {input_file_base_name}: {str(e)}")
        return False  # Failure

# List of gRNA sequences to highlight
grna_list = ["CAGCAGCTGG","TCCAACCAGG"]

# Initialize counters and timer
file_count = 0
successful_files = 0
start_time = time.time()
start_timestamp = datetime.fromtimestamp(start_time).strftime('%H:%M:%S')

# Process all .ab1 files
for file_name in os.listdir(input_folder_path):
    if file_name.endswith(".ab1"):
        file_count += 1
        input_file_path = os.path.join(input_folder_path, file_name)
        try:
            success = process_input_file(input_file_path, grna_sequences=grna_list, driver=driver)
            if not success:
                # Reinitialize driver on failure
                driver.quit()
                driver = init_driver()
                # Retry the file
                success = process_input_file(input_file_path, grna_sequences=grna_list, driver=driver)
            if success:
                successful_files += 1
        except InvalidSessionIdException:
            print(f"Session crashed for {file_name}. Reinitializing driver.")
            driver.quit()
            driver = init_driver()
            success = process_input_file(input_file_path, grna_sequences=grna_list, driver=driver)
            if success:
                successful_files += 1

# Calculate total time taken and end timestamp
end_time = time.time()
end_timestamp = datetime.fromtimestamp(end_time).strftime('%H:%M:%S')
total_time = end_time - start_time

# Print and log summary
summary = (
    f"\nAnalysis Summary:\n"
    f"Total files analyzed: {file_count}\n"
    f"Successfully processed files: {successful_files}\n"
    f"Start time: {start_timestamp}\n"
    f"End time: {end_timestamp}\n"
    f"Total time taken: {total_time:.2f} seconds\n"
)
print(summary)
with open(log_file_path, "a", encoding="utf-8") as log_file:
    log_file.write(summary)

driver.quit()

© 2025 Tata Institute for Genetics and Society (TIGS), Bangalore, India
All rights reserved.
