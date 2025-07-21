import json
import time
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.chrome.options import Options
from selenium.common.exceptions import TimeoutException, NoSuchElementException
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class AbortionFinderScraper:
    def __init__(self, headless=True):
        """Initialize the scraper with Chrome options"""
        chrome_options = Options()
        if headless:
            chrome_options.add_argument("--headless")
        chrome_options.add_argument("--no-sandbox")
        chrome_options.add_argument("--disable-dev-shm-usage")
        chrome_options.add_argument("--user-agent=Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36")
        
        self.driver = webdriver.Chrome(options=chrome_options)
        self.wait = WebDriverWait(self.driver, 10)
        self.providers = []
    
    def scrape_provider_info(self, provider_element):
        """Extract information from a single provider element"""
        try:
            provider_info = {}
            
            # Extract name
            try:
                name_element = provider_element.find_element(By.CSS_SELECTOR, "h3, .provider-name, [class*='name']")
                provider_info['name'] = name_element.text.strip()
            except NoSuchElementException:
                provider_info['name'] = "Name not found"
            
            # Extract address
            try:
                address_element = provider_element.find_element(By.CSS_SELECTOR, ".address, [class*='address']")
                provider_info['address'] = address_element.text.strip()
            except NoSuchElementException:
                provider_info['address'] = "Address not found"
            
            # Extract phone number
            try:
                phone_element = provider_element.find_element(By.CSS_SELECTOR, ".phone, [class*='phone'], a[href^='tel:']")
                provider_info['phone'] = phone_element.text.strip()
            except NoSuchElementException:
                provider_info['phone'] = "Phone not found"
            
            # Extract website/link
            try:
                link_element = provider_element.find_element(By.CSS_SELECTOR, "a[href^='http'], .website, [class*='website']")
                provider_info['website'] = link_element.get_attribute('href')
            except NoSuchElementException:
                provider_info['website'] = "Website not found"
            
            # Extract any additional information (services, hours, etc.)
            try:
                info_elements = provider_element.find_elements(By.CSS_SELECTOR, ".info, .details, [class*='info'], [class*='detail']")
                additional_info = []
                for element in info_elements:
                    text = element.text.strip()
                    if text and text not in [provider_info.get('name', ''), provider_info.get('address', ''), provider_info.get('phone', '')]:
                        additional_info.append(text)
                provider_info['additional_info'] = additional_info
            except NoSuchElementException:
                provider_info['additional_info'] = []
            
            return provider_info
            
        except Exception as e:
            logger.error(f"Error extracting provider info: {str(e)}")
            return None
    
    def scrape_current_page(self):
        """Scrape all providers from the current page"""
        try:
            # Wait for providers to load
            self.wait.until(EC.presence_of_element_located((By.CSS_SELECTOR, "[class*='provider'], .clinic, .facility")))
            
            # Find all provider elements - try different selectors
            provider_selectors = [
                "[class*='provider']",
                ".clinic",
                ".facility", 
                "[class*='clinic']",
                "[class*='facility']",
                ".result-item",
                "[class*='result']"
            ]
            
            providers_found = []
            for selector in provider_selectors:
                try:
                    elements = self.driver.find_elements(By.CSS_SELECTOR, selector)
                    if elements:
                        providers_found = elements
                        logger.info(f"Found {len(elements)} providers using selector: {selector}")
                        break
                except:
                    continue
            
            if not providers_found:
                # Fallback: try to find any clickable elements or cards
                providers_found = self.driver.find_elements(By.CSS_SELECTOR, "[class*='card'], .item, [role='article']")
            
            page_providers = []
            for provider_element in providers_found:
                provider_info = self.scrape_provider_info(provider_element)
                if provider_info:
                    page_providers.append(provider_info)
            
            logger.info(f"Scraped {len(page_providers)} providers from current page")
            return page_providers
            
        except TimeoutException:
            logger.error("Timeout waiting for providers to load")
            return []
        except Exception as e:
            logger.error(f"Error scraping current page: {str(e)}")
            return []
    
    def get_next_page_button(self):
        """Find and return the next page button"""
        next_button_selectors = [
            "a[aria-label='Next']",
            ".next",
            "[class*='next']",
            ".pagination .next",
            "a[rel='next']",
            ".page-numbers .next"
        ]
        
        for selector in next_button_selectors:
            try:
                button = self.driver.find_element(By.CSS_SELECTOR, selector)
                if button.is_enabled() and button.is_displayed():
                    return button
            except NoSuchElementException:
                continue
        
        # Try to find pagination numbers and get the next one
        try:
            current_page = self.driver.find_element(By.CSS_SELECTOR, ".current, .active, [aria-current='page']")
            current_num = int(current_page.text.strip())
            next_page_link = self.driver.find_element(By.XPATH, f"//a[text()='{current_num + 1}']")
            return next_page_link
        except:
            pass
        
        return None
    
    def scrape_all_pages(self, base_url, max_pages=None):
        """Scrape all pages of providers"""
        try:
            self.driver.get(base_url)
            logger.info(f"Starting scrape of: {base_url}")
            
            # Wait for initial page load
            time.sleep(3)
            
            page_num = 1
            all_providers = []
            
            while True:
                logger.info(f"Scraping page {page_num}")
                
                # Scrape current page
                page_providers = self.scrape_current_page()
                all_providers.extend(page_providers)
                
                # Check if we've reached max pages
                if max_pages and page_num >= max_pages:
                    logger.info(f"Reached maximum pages limit: {max_pages}")
                    break
                
                # Look for next page button
                next_button = self.get_next_page_button()
                
                if next_button:
                    try:
                        # Scroll to button and click
                        self.driver.execute_script("arguments[0].scrollIntoView();", next_button)
                        time.sleep(1)
                        next_button.click()
                        
                        # Wait for new page to load
                        time.sleep(3)
                        page_num += 1
                        
                    except Exception as e:
                        logger.error(f"Error clicking next page: {str(e)}")
                        break
                else:
                    logger.info("No more pages found")
                    break
            
            self.providers = all_providers
            logger.info(f"Total providers scraped: {len(all_providers)}")
            return all_providers
            
        except Exception as e:
            logger.error(f"Error during scraping: {str(e)}")
            return []
    
    def save_to_json(self, filename="abortion_providers_california.json"):
        """Save scraped providers to JSON file"""
        try:
            with open(filename, 'w', encoding='utf-8') as f:
                json.dump(self.providers, f, indent=2, ensure_ascii=False)
            logger.info(f"Data saved to {filename}")
        except Exception as e:
            logger.error(f"Error saving to JSON: {str(e)}")
    
    def close(self):
        """Close the browser"""
        if self.driver:
            self.driver.quit()

def main():
    """Main function to run the scraper"""
    url = "https://www.abortionfinder.org/abortion-guides-by-state/abortion-in-california/providers"
    
    scraper = AbortionFinderScraper(headless=False)  # Set to True for headless mode
    
    try:
        # Scrape all providers
        providers = scraper.scrape_all_pages(url, max_pages=19)  # Limit to 10 pages, remove limit by setting to None
        
        # Save to JSON
        if providers:
            scraper.save_to_json("california_abortion_providers.json")
            
            # Print summary
            print(f"\nScraping completed!")
            print(f"Total providers found: {len(providers)}")
            if providers:
                print("\nSample provider:")
                print(json.dumps(providers[0], indent=2))
        else:
            print("No providers were found. The website structure may have changed.")
            
    except Exception as e:
        logger.error(f"Script error: {str(e)}")
    finally:
        scraper.close()

if __name__ == "__main__":
    main()