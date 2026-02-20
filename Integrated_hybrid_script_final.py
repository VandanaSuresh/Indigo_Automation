# Automated Sanger sequence Data Analysis INDIGO & ICE 
#Created by Vandana Suresh
# Research assistant @ TIGS, Bangalore
# Date-05/02/2026

import os
import time
import re
import sys
import json
import pandas as pd
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from datetime import datetime
from selenium.common.exceptions import (
    InvalidSessionIdException, TimeoutException, NoSuchElementException,
    WebDriverException, StaleElementReferenceException
)
from Bio.Seq import Seq
from Bio import SeqIO

# ==================================================================================
# CUSTOM ERROR CLASSES
# ==================================================================================
class AnalysisError(Exception):
    """Base exception for analysis errors"""
    pass

class IndigoError(AnalysisError):
    """INDIGO-specific errors"""
    pass

class ICEError(AnalysisError):
    """ICE-specific errors"""
    pass

class PrerequisiteError(AnalysisError):
    """Prerequisites validation errors"""
    pass

# ==================================================================================
# LOGGER CLASS - BETTER LOGGING
# ==================================================================================
class AnalysisLogger:
    """Centralized logging for console and file"""
    
    def __init__(self, log_file_path):
        self.log_file_path = log_file_path
        self.errors = []
        self.warnings = []
        
    def info(self, message):
        """Log info message"""
        print(message)
        self._write_to_file(f"[INFO] {message}")
    
    def error(self, message):
        """Log error message"""
        print(f"✗ ERROR: {message}")
        self._write_to_file(f"[ERROR] {message}")
        self.errors.append(message)
    
    def warning(self, message):
        """Log warning message"""
        print(f"⚠ WARNING: {message}")
        self._write_to_file(f"[WARNING] {message}")
        self.warnings.append(message)
    
    def success(self, message):
        """Log success message"""
        print(f"✓ {message}")
        self._write_to_file(f"[SUCCESS] {message}")
    
    def debug(self, message):
        """Log debug message"""
        self._write_to_file(f"[DEBUG] {message}")
    
    def _write_to_file(self, message):
        """Write message to log file with error handling"""
        try:
            with open(self.log_file_path, 'a', encoding='utf-8') as f:
                timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
                f.write(f"{timestamp} - {message}\n")
        except IOError as e:
            print(f"⚠ Could not write to log file: {e}")
        except Exception as e:
            print(f"⚠ Unexpected logging error: {e}")
    
    def get_summary(self):
        """Get error and warning summary"""
        return {
            'error_count': len(self.errors),
            'warning_count': len(self.warnings),
            'errors': self.errors[:10],  # Last 10 errors
            'warnings': self.warnings[:10]  # Last 10 warnings
        }

# ==================================================================================
# BIOPYTHON COMPATIBILITY PATCH FOR ICE
# ==================================================================================
try:
    from Bio.Align import MultipleSeqAlignment
    from Bio import AlignIO
    from io import StringIO

    def alignment_format_patch(self, format):
        handle = StringIO()
        AlignIO.write(self, handle, format)
        return handle.getvalue()

    MultipleSeqAlignment.format = alignment_format_patch
    print("✓ Biopython compatibility patch applied.")
except Exception as e:
    print(f"⚠ Warning: Could not apply Biopython patch: {e}")

# ==================================================================================
# CONFIGURATION - UPDATE THESE PATHS
# ==================================================================================
# Input/Output Directories
input_folder_path = r""
wild_type_file_path = r""
download_dir = r""
log_file_path = os.path.join(download_dir, "analysis_log.txt")

# ChromeDriver path
chromedriver_path = r"C:\New_folder\chromedriver-win64\chromedriver-win64\chromedriver.exe"

# ICE module path
ice_source_path = r"C:\Users\Vandana.S\Desktop\TIGS\Automation_new\HYBRID_WORK\ice-master"

# Analysis Parameters
grna_sequences = [""]
ICE_TARGET_SEQUENCE_FALLBACK = ""

# Retry Configuration
MAX_RETRIES_INDIGO = 2  # Max retries for INDIGO analysis
MAX_RETRIES_ICE = 1  # Max retries for ICE analysis
INDIGO_WAIT_TIME = 15  # Seconds to wait for INDIGO analysis
DRIVER_TIMEOUT = 30  # Selenium WebDriver timeout

# ==================================================================================
# CREATE OUTPUT DIRECTORY STRUCTURE
# ==================================================================================
try:
    os.makedirs(download_dir, exist_ok=True)
    indigo_output_dir = download_dir
    ice_output_dir = os.path.join(os.path.dirname(download_dir), "ICE_RESULTS")
    os.makedirs(ice_output_dir, exist_ok=True)
except OSError as e:
    print(f"✗ Failed to create output directories: {e}")
    sys.exit(1)

# Initialize logger
logger = AnalysisLogger(log_file_path)
logger.info(f"Log file created at: {log_file_path}")

# ==================================================================================
# VALIDATION & SETUP
# ==================================================================================
def validate_prerequisites():
    """Validate all prerequisites before starting with detailed error handling"""
    logger.info("\n" + "="*80)
    logger.info("VALIDATING PREREQUISITES")
    logger.info("="*80)
    
    try:
        # Check input folder
        if not os.path.exists(input_folder_path):
            raise PrerequisiteError(f"Input folder not found: {input_folder_path}")
        
        if not os.path.isdir(input_folder_path):
            raise PrerequisiteError(f"Input path is not a directory: {input_folder_path}")
        
        logger.success(f"Input folder exists: {input_folder_path}")
        
        # Check wildtype file
        if not os.path.exists(wild_type_file_path):
            raise PrerequisiteError(f"Wildtype file not found: {wild_type_file_path}")
        
        if not os.path.isfile(wild_type_file_path):
            raise PrerequisiteError(f"Wildtype path is not a file: {wild_type_file_path}")
        
        # Try to read wildtype file
        try:
            record = SeqIO.read(wild_type_file_path, "abi")
            logger.success(f"Wildtype file valid (.ab1 format): {os.path.basename(wild_type_file_path)}")
        except Exception as e:
            raise PrerequisiteError(f"Wildtype file cannot be read as .ab1: {e}")
        
        # Check ChromeDriver
        if not os.path.exists(chromedriver_path):
            raise PrerequisiteError(f"ChromeDriver not found: {chromedriver_path}")
        
        if not os.path.isfile(chromedriver_path):
            raise PrerequisiteError(f"ChromeDriver path is not a file: {chromedriver_path}")
        
        logger.success(f"ChromeDriver found")
        
        # Check input files
        try:
            all_files = os.listdir(input_folder_path)
            ab1_files = [f for f in all_files if f.endswith('.ab1')]
        except OSError as e:
            raise PrerequisiteError(f"Cannot read input folder: {e}")
        
        if not ab1_files:
            raise PrerequisiteError(f"No .ab1 files found in: {input_folder_path}")
        
        logger.success(f"Found {len(ab1_files)} .ab1 files to process")
        
        # Verify each file can be read
        unreadable_files = []
        for ab1_file in ab1_files:
            file_path = os.path.join(input_folder_path, ab1_file)
            try:
                SeqIO.read(file_path, "abi")
            except Exception as e:
                unreadable_files.append((ab1_file, str(e)))
        
        if unreadable_files:
            logger.warning(f"{len(unreadable_files)} file(s) cannot be read as .ab1:")
            for fname, error in unreadable_files[:5]:
                logger.warning(f"  - {fname}: {error[:50]}")
            # Continue anyway, but track these
        
        logger.info("="*80 + "\n")
        return ab1_files
        
    except PrerequisiteError as e:
        logger.error(f"Prerequisite validation failed: {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Unexpected error during prerequisite validation: {e}")
        import traceback
        logger.debug(traceback.format_exc())
        sys.exit(1)

# ==================================================================================
# LOAD ICE MODULE
# ==================================================================================
def load_ice_module():
    """Load ICE module with error handling"""
    global ICE_AVAILABLE
    
    if ice_source_path not in sys.path:
        sys.path.insert(0, ice_source_path)
    
    try:
        from ice.analysis import single_sanger_analysis
        logger.success("ICE Module loaded successfully")
        ICE_AVAILABLE = True
        return single_sanger_analysis
    except ImportError as e:
        logger.warning(f"ICE Module import failed: {e}")
        logger.warning("Analysis will continue with INDIGO only (no ICE fallback)")
        ICE_AVAILABLE = False
        return None
    except Exception as e:
        logger.warning(f"Unexpected error loading ICE module: {e}")
        ICE_AVAILABLE = False
        return None

single_sanger_analysis = load_ice_module()
ICE_AVAILABLE = single_sanger_analysis is not None

# ==================================================================================
# SELENIUM FUNCTIONS - IMPROVED ERROR HANDLING
# ==================================================================================
def init_driver():
    """Initialize WebDriver with error handling"""
    try:
        service = Service(chromedriver_path)
        service.log_path = os.devnull
        options = webdriver.ChromeOptions()
        
        prefs = {"download.default_directory": indigo_output_dir}
        options.add_experimental_option("prefs", prefs)
        options.add_argument("--log-level=3")
        options.add_experimental_option("excludeSwitches", ["enable-logging"])
        # options.add_argument("--headless=new")
        
        driver = webdriver.Chrome(service=service, options=options)
        driver.set_page_load_timeout(DRIVER_TIMEOUT)
        driver.set_script_timeout(DRIVER_TIMEOUT)
        
        return driver
        
    except WebDriverException as e:
        raise IndigoError(f"Failed to initialize WebDriver: {e}")
    except Exception as e:
        raise IndigoError(f"Unexpected error initializing WebDriver: {e}")

def highlight_pam_sequence(html_content, grna_sequences=None):
    """Highlight gRNA sequences with error handling"""
    try:
        if grna_sequences is None:
            grna_sequences = [""]
        
        if not grna_sequences:
            logger.warning("No gRNA sequences provided for highlighting")
            return html_content
        
        colors = ["cyan"]
        for i, grna_sequence in enumerate(grna_sequences):
            if not grna_sequence:
                logger.warning(f"Empty gRNA sequence at index {i}")
                continue
            
            color = colors[i % len(colors)]
            pattern = re.compile(re.escape(grna_sequence), re.IGNORECASE)
            html_content = pattern.sub(
                lambda m: f'<span style="background-color: {color}; font-weight: bold;">{m.group(0)}</span>',
                html_content
            )
        
        return html_content
        
    except Exception as e:
        logger.warning(f"Error highlighting PAM sequence: {e}")
        return html_content  # Return original if highlighting fails

def process_input_file(input_file_path, grna_sequences, driver, retry_count=0):
    """
    Process file with INDIGO with comprehensive error handling
    """
    input_file_base_name = os.path.splitext(os.path.basename(input_file_path))[0]
    
    try:
        # Verify input file exists
        if not os.path.exists(input_file_path):
            raise IndigoError(f"Input file not found: {input_file_path}")
        
        # Load INDIGO page
        try:
            driver.get("https://www.gear-genomics.com/indigo/")
        except TimeoutException:
            raise IndigoError("Timeout loading INDIGO website")
        except WebDriverException as e:
            raise IndigoError(f"WebDriver error accessing INDIGO: {e}")
        
        wait = WebDriverWait(driver, DRIVER_TIMEOUT)
        
        # Upload input file
        try:
            input_file_upload = wait.until(EC.presence_of_element_located((By.ID, "inputFile")))
            input_file_upload.send_keys(input_file_path)
            logger.debug(f"Uploaded {input_file_base_name}")
        except TimeoutException:
            raise IndigoError("Timeout waiting for input file upload element")
        except StaleElementReferenceException:
            raise IndigoError("Element reference became stale during upload")
        except Exception as e:
            raise IndigoError(f"Error uploading input file: {e}")
        
        # Verify upload
        time.sleep(2)
        try:
            file_value = input_file_upload.get_attribute("value")
            if not file_value:
                raise IndigoError("Sample file upload verification failed - no file selected")
        except Exception as e:
            raise IndigoError(f"Error verifying file upload: {e}")
        
        # Click wildtype tab
        try:
            wild_type_link = wait.until(EC.element_to_be_clickable((By.ID, "target-chromatogram-tab")))
            driver.execute_script("arguments[0].scrollIntoView();", wild_type_link)
            time.sleep(1)
            driver.execute_script("arguments[0].click();", wild_type_link)
        except TimeoutException:
            raise IndigoError("Timeout clicking wildtype chromatogram tab")
        except Exception as e:
            raise IndigoError(f"Error clicking wildtype tab: {e}")
        
        # Upload wildtype file
        try:
            wild_type_file_input = wait.until(EC.presence_of_element_located((By.ID, "targetFileChromatogram")))
            wild_type_file_input.send_keys(wild_type_file_path)
            logger.debug(f"Uploaded wildtype file")
        except TimeoutException:
            raise IndigoError("Timeout waiting for wildtype file upload element")
        except Exception as e:
            raise IndigoError(f"Error uploading wildtype file: {e}")
        
        # Verify wildtype upload
        time.sleep(2)
        try:
            wt_value = wild_type_file_input.get_attribute("value")
            if not wt_value:
                raise IndigoError("Wildtype file upload verification failed")
        except Exception as e:
            raise IndigoError(f"Error verifying wildtype upload: {e}")
        
        # Submit analysis
        try:
            submit_button = wait.until(EC.element_to_be_clickable((By.ID, "btn-submit")))
            driver.execute_script("arguments[0].scrollIntoView();", submit_button)
            time.sleep(1)
            driver.execute_script("arguments[0].click();", submit_button)
        except TimeoutException:
            raise IndigoError("Timeout clicking submit button")
        except Exception as e:
            raise IndigoError(f"Error submitting analysis: {e}")
        
        # Wait for INDIGO to process
        logger.debug(f"Waiting {INDIGO_WAIT_TIME} seconds for INDIGO analysis...")
        time.sleep(INDIGO_WAIT_TIME)
        
        # Try to get results
        try:
            result_download_link = wait.until(
                EC.presence_of_element_located((By.LINK_TEXT, "Download HTML")),
                message="Download link not found"
            )
            driver.execute_script("arguments[0].click();", result_download_link)
            time.sleep(5)
            logger.debug(f"Downloaded HTML for {input_file_base_name}")
            return True, None
            
        except (TimeoutException, NoSuchElementException):
            # Fallback to page source
            try:
                page_source = driver.page_source
                
                # Check for errors
                error_patterns = [
                    ("Error in running Indigo:", r'Error in running Indigo: ([^<]+)'),
                    ("Alignment of trace to reference failed", None),
                    ("execution halted", None),
                    ("package:stats", None),
                ]
                
                for error_text, error_pattern in error_patterns:
                    if error_text in page_source:
                        if error_pattern:
                            match = re.search(error_pattern, page_source)
                            if match:
                                return False, match.group(1)
                        return False, error_text
                
                # No errors detected, save page source
                highlighted_html = highlight_pam_sequence(page_source, grna_sequences=grna_sequences)
                result_html_file_path = os.path.join(indigo_output_dir, f"{input_file_base_name}_results_highlighted.html")
                
                with open(result_html_file_path, "w", encoding="utf-8") as html_file:
                    html_file.write(highlighted_html)
                
                logger.debug(f"Saved results (page source) for {input_file_base_name}")
                return True, None
                
            except IOError as e:
                raise IndigoError(f"Error saving results to file: {e}")
            except Exception as e:
                raise IndigoError(f"Error processing results: {e}")
        
    except IndigoError as e:
        logger.debug(f"INDIGO Error: {str(e)}")
        return False, str(e)
    except StaleElementReferenceException:
        if retry_count < MAX_RETRIES_INDIGO:
            logger.warning(f"Stale element reference, retrying (attempt {retry_count + 1})")
            time.sleep(2)
            return process_input_file(input_file_path, grna_sequences, driver, retry_count + 1)
        else:
            return False, "Stale element reference after retries"
    except WebDriverException as e:
        return False, f"WebDriver error: {str(e)[:100]}"
    except Exception as e:
        logger.debug(f"Unexpected error in process_input_file: {e}")
        import traceback
        logger.debug(traceback.format_exc())
        return False, f"Unexpected error: {str(e)[:100]}"

def process_with_ice(input_file_path, ice_target_sequence, retry_count=0):
    """Process file with ICE as fallback with error handling"""
    
    if not ICE_AVAILABLE or single_sanger_analysis is None:
        return False, None, "ICE module not available"
    
    file_name = os.path.basename(input_file_path)
    sample_name = os.path.splitext(file_name)[0]
    sample_dir = os.path.join(ice_output_dir, sample_name)
    
    try:
        # Create sample directory
        try:
            os.makedirs(sample_dir, exist_ok=True)
        except OSError as e:
            raise ICEError(f"Cannot create output directory: {e}")
        
        # Verify input file
        if not os.path.exists(input_file_path):
            raise ICEError(f"Input file not found: {input_file_path}")
        
        if not os.path.exists(wild_type_file_path):
            raise ICEError(f"Wildtype file not found: {wild_type_file_path}")
        
        # Verify target sequence
        if not ice_target_sequence or not isinstance(ice_target_sequence, str):
            raise ICEError(f"Invalid target sequence: {ice_target_sequence}")
        
        # Run ICE analysis
        try:
            result_json = single_sanger_analysis(
                control_path=wild_type_file_path,
                sample_path=input_file_path,
                base_outputname=os.path.join(sample_dir, "ICE"),
                guide=ice_target_sequence,
                verbose=False
            )
        except Exception as e:
            error_msg = str(e)
            if "not found in control sequence" in error_msg:
                raise ICEError(f"Target sequence not found in control sequence")
            elif "No such file or directory" in error_msg:
                raise ICEError(f"File not found during ICE analysis")
            else:
                raise ICEError(f"ICE analysis failed: {error_msg[:100]}")
        
        # Parse results
        try:
            if isinstance(result_json, str):
                res_data = json.loads(result_json)
            else:
                res_data = result_json
            
            indel = res_data.get('ice', 0)
            r2 = res_data.get('rsq', 0)
            
            # Validate results
            if indel is None:
                indel = 0.0
            if r2 is None:
                r2 = 0.0
            
            indel = float(indel)
            r2 = float(r2)
            
            # Check for suspicious results
            if indel < 0 or indel > 100:
                logger.warning(f"Unusual indel percentage for {sample_name}: {indel}%")
            
            if r2 < 0 or r2 > 1:
                logger.warning(f"Unusual R² value for {sample_name}: {r2}")
            
            result_dict = {
                'indel_percentage': indel,
                'r_squared': r2
            }
            
            return True, result_dict, None
            
        except (json.JSONDecodeError, ValueError, KeyError) as e:
            raise ICEError(f"Error parsing ICE results: {e}")
        
    except ICEError as e:
        logger.debug(f"ICE Error: {str(e)}")
        return False, None, str(e)
    except Exception as e:
        logger.debug(f"Unexpected error in process_with_ice: {e}")
        import traceback
        logger.debug(traceback.format_exc())
        return False, None, f"Unexpected error: {str(e)[:100]}"

# ==================================================================================
# MAIN ANALYSIS WORKFLOW
# ==================================================================================
def main():
    """Main analysis workflow with comprehensive error handling"""
    
    # Validate prerequisites
    ab1_files = validate_prerequisites()
    
    # Initialize tracking
    file_count = 0
    successful_files = 0
    failed_files = 0
    results_tracker = []
    start_time = time.time()
    start_timestamp = datetime.fromtimestamp(start_time).strftime('%Y-%m-%d %H:%M:%S')
    
    total_files = len(ab1_files)
    
    logger.info("="*80)
    logger.info("SANGER SEQUENCING HYBRID ANALYSIS (IMPROVED ERROR HANDLING)")
    logger.info("="*80)
    logger.info(f"Total files to process: {total_files}")
    logger.info(f"Start time: {start_timestamp}")
    logger.info(f"gRNA sequences: {grna_sequences}")
    logger.info(f"ICE Available: {ICE_AVAILABLE}")
    logger.info("="*80 + "\n")
    
    # Initialize driver
    driver = None
    try:
        driver = init_driver()
        logger.success("WebDriver initialized successfully")
    except IndigoError as e:
        logger.error(f"Failed to initialize WebDriver: {e}")
        logger.info("Cannot continue without WebDriver. Exiting.")
        sys.exit(1)
    
    # Process each file
    for file_name in sorted(ab1_files):
        file_count += 1
        input_file_path = os.path.join(input_folder_path, file_name)
        
        logger.info(f"[{file_count}/{total_files}] Processing: {file_name}")
        
        try:
            # Try INDIGO
            success, indigo_error = process_input_file(
                input_file_path,
                grna_sequences=grna_sequences,
                driver=driver
            )
            
            if success:
                logger.success(f"INDIGO analysis successful")
                successful_files += 1
                results_tracker.append({
                    'Sample': file_name,
                    'Primary_Tool': 'Indigo',
                    'Status': 'Success',
                    'Fallback_Used': 'No',
                    'Error': ''
                })
                continue
            else:
                logger.warning(f"INDIGO analysis failed: {indigo_error[:60] if indigo_error else 'Unknown'}")
            
        except InvalidSessionIdException as e:
            logger.warning("WebDriver session crashed. Reinitializing...")
            try:
                if driver:
                    driver.quit()
                driver = init_driver()
                
                # Retry INDIGO
                success, indigo_error = process_input_file(
                    input_file_path,
                    grna_sequences=grna_sequences,
                    driver=driver
                )
                
                if success:
                    logger.success(f"INDIGO analysis successful (after retry)")
                    successful_files += 1
                    results_tracker.append({
                        'Sample': file_name,
                        'Primary_Tool': 'Indigo',
                        'Status': 'Success',
                        'Fallback_Used': 'No',
                        'Error': ''
                    })
                    continue
                    
            except Exception as retry_error:
                logger.error(f"Failed to reinitialize driver: {retry_error}")
                indigo_error = str(retry_error)
        
        except Exception as e:
            logger.error(f"Unexpected error during INDIGO processing: {e}")
            indigo_error = str(e)
        
        # Route to ICE fallback
        if ICE_AVAILABLE:
            logger.info("Routing to ICE fallback...")
            
            try:
                ice_success, ice_results, ice_error = process_with_ice(
                    input_file_path,
                    ICE_TARGET_SEQUENCE_FALLBACK
                )
                
                if ice_success and ice_results:
                    logger.success(f"ICE analysis successful")
                    logger.info(f"  Indel: {ice_results['indel_percentage']:.2f}% | R²: {ice_results['r_squared']:.4f}")
                    successful_files += 1
                    results_tracker.append({
                        'Sample': file_name,
                        'Primary_Tool': 'Indigo',
                        'Status': 'Success (via ICE)',
                        'Fallback_Used': 'Yes',
                        'Indigo_Error': indigo_error[:80] if indigo_error else '',
                        'ICE_Indel_%': f"{ice_results['indel_percentage']:.2f}",
                        'ICE_R²': f"{ice_results['r_squared']:.4f}",
                        'Error': ''
                    })
                else:
                    logger.error(f"ICE analysis failed: {ice_error[:60] if ice_error else 'Unknown'}")
                    failed_files += 1
                    results_tracker.append({
                        'Sample': file_name,
                        'Primary_Tool': 'Indigo',
                        'Status': 'Failed (both tools)',
                        'Fallback_Used': 'Attempted',
                        'Indigo_Error': indigo_error[:80] if indigo_error else '',
                        'ICE_Error': ice_error[:80] if ice_error else '',
                        'Error': ice_error[:80] if ice_error else ''
                    })
                
            except Exception as e:
                logger.error(f"Unexpected error in ICE fallback: {e}")
                failed_files += 1
                results_tracker.append({
                    'Sample': file_name,
                    'Primary_Tool': 'Indigo',
                    'Status': 'Failed (both tools)',
                    'Fallback_Used': 'Attempted',
                    'Indigo_Error': indigo_error[:80] if indigo_error else '',
                    'ICE_Error': str(e)[:80],
                    'Error': str(e)[:80]
                })
        
        else:
            logger.warning("ICE not available, skipping fallback")
            failed_files += 1
            results_tracker.append({
                'Sample': file_name,
                'Primary_Tool': 'Indigo',
                'Status': 'Failed (ICE unavailable)',
                'Fallback_Used': 'N/A',
                'Indigo_Error': indigo_error[:80] if indigo_error else '',
                'Error': 'ICE module not available'
            })
        
        logger.info("")  # Blank line for readability
    
    # Cleanup
    if driver:
        try:
            driver.quit()
            logger.success("WebDriver closed successfully")
        except Exception as e:
            logger.warning(f"Error closing WebDriver: {e}")
    
    # Generate report
    try:
        end_time = time.time()
        end_timestamp = datetime.fromtimestamp(end_time).strftime('%Y-%m-%d %H:%M:%S')
        total_time = end_time - start_time
        
        # Create DataFrame
        df = pd.DataFrame(results_tracker)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        report_path = os.path.join(os.path.dirname(indigo_output_dir), f"hybrid_analysis_report_{timestamp}.csv")
        
        df.to_csv(report_path, index=False)
        logger.success(f"Report saved: {report_path}")
        
    except Exception as e:
        logger.error(f"Error generating report: {e}")
        report_path = "Unable to save"
    
    # Print summary
    summary = (
        f"\n{'='*80}\n"
        f"ANALYSIS COMPLETE\n"
        f"{'='*80}\n"
        f"Total files analyzed: {file_count}\n"
        f"Successfully processed: {successful_files}\n"
        f"Failed: {failed_files}\n"
        f"Success rate: {successful_files/file_count*100:.1f}%\n"
        f"Start time: {start_timestamp}\n"
        f"End time: {end_timestamp}\n"
        f"Total time: {total_time:.2f} seconds ({total_time/60:.1f} minutes)\n"
        f"Report: {report_path}\n"
    )
    
    # Add error summary if there were issues
    error_summary = logger.get_summary()
    if error_summary['error_count'] > 0 or error_summary['warning_count'] > 0:
        summary += (
            f"\nIssues Encountered:\n"
            f"  Errors: {error_summary['error_count']}\n"
            f"  Warnings: {error_summary['warning_count']}\n"
        )
    
    summary += f"{'='*80}\n"
    
    logger.info(summary)

# ==================================================================================
# RUN SCRIPT
# ==================================================================================
if __name__ == "__main__":
    try:
        main()
        logger.success("PIPELINE COMPLETED SUCCESSFULLY")
    except KeyboardInterrupt:
        logger.warning("Pipeline interrupted by user")
        sys.exit(0)
    except Exception as e:
        logger.error(f"CRITICAL PIPELINE ERROR: {e}")
        import traceback
        logger.debug(traceback.format_exc())
        sys.exit(1)