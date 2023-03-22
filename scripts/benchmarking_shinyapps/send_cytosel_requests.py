#!/usr/bin/env python

import time
from selenium.webdriver.common.by import By
from selenium import webdriver
from selenium.webdriver.chrome.options import Options as ChromeOptions
from selenium.webdriver.support.wait import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import TimeoutException, ElementClickInterceptedException, NoAlertPresentException

chrome_options = ChromeOptions()
chrome_options.add_argument('--headless')
chrome_options.add_argument('--disable-gpu')

driver = webdriver.Chrome(options=chrome_options)
iterations = 0

while iterations < 10:
    try:
        driver.get('https://camlab.shinyapps.io/cytosel/')
        try:
            alert = driver.switch_to.alert
            alert.accept()
        except NoAlertPresentException:
            pass

        WebDriverWait(driver, 180).until(
            EC.visibility_of_element_located((By.ID, "curated_dataset"))
        )

        try:
            alert = driver.switch_to.alert
            alert.accept()
        except NoAlertPresentException:
            pass

        driver.find_element(By.ID, "curated_dataset").click()
        WebDriverWait(driver, 180).until(
            EC.visibility_of_element_located((By.ID, "pick_curated"))
        )
        driver.find_element(By.ID, "pick_curated").click()
        print("curated selected")

        element = WebDriverWait(driver, 180).until(
            EC.visibility_of_element_located((By.ID, "start_analysis")))

        element.click()

        try:
            alert = driver.switch_to.alert
            alert.accept()
        except NoAlertPresentException:
            pass

        print("starting analysis")

        # gene_tab_visible = WebDriverWait(driver, 180).until(
        #         EC.visibility_of_element_located((By.XPATH, '//a[contains(@href,"#shiny-tab-gene-expression")]'))
        #     )

        analysis_successful = WebDriverWait(driver, 300).until(
            EC.visibility_of_element_located((By.ID, "markers_change_modal")))

        print("panel generated.")

        umap_is_ready = WebDriverWait(driver, 180).until(
            EC.visibility_of_element_located((By.XPATH, '//a[contains(@href,"#shiny-tab-UMAP")]')))
        umap_is_ready.click()

        # umap_visible = WebDriverWait(driver, 180).until(
        #         EC.visibility_of_element_located((By.ID, "show_umap_legend")))

        umap_is_visible = WebDriverWait(driver, 180).until(
            EC.visibility_of_element_located((By.XPATH, "//span[.='Show UMAP plot legends']")))

        time.sleep(3)

        print("rendering UMAP.")

    except (TimeoutException, ElementClickInterceptedException) as error:
        print("Unable to perform the complete cytosel analysis.")

    iterations += 1



