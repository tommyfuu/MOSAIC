__author__ = "Tom Fu"

__version__ = "1.0"

# next steps: take in more different formats
# (1) excel, txt, etc.
# (2) take in non empty matrix, instead two lists, should be straightforward

import pandas as pd
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
import os
import time
from webdriver_manager.chrome import ChromeDriverManager
from selenium.webdriver.common.action_chains import ActionChains
from selenium.common.exceptions import TimeoutException
from selenium.webdriver.common.by import By


def pubtatorAutoSearch(meshTerm, headless):
    """fetch relevant stats for one particular meshterm and return the search results"""
    options = webdriver.chrome.options.Options()
    # get window of the right size
    options.add_argument("window-size=1200x600")
    options.headless = headless
    driver = webdriver.Chrome()
    driver.set_page_load_timeout(8)

    # get website
    # handles the long pending web page loading
    loaded = False
    while loaded == False:
        try:
            driver.get('https://www.ncbi.nlm.nih.gov/research/pubtator3/')
            loaded = True
        except TimeoutException:
            driver.close()
            driver = webdriver.Chrome(
                ChromeDriverManager().install(), options=options)
            driver.set_page_load_timeout(8)

    search_input = driver.find_element(By.CSS_SELECTOR,
        'body > app-root > div > app-home > div.container.litsearch-intro > div > div.intro-main > div > div.col.l8.m10.offset-m1.s12 > div:nth-child(2) > app-search > div > div > div.input-field > div > div > input')
    print('aaa', search_input)
    time.sleep(1)

    search_input.send_keys(meshTerm+"\n")
    time.sleep(8)

    # get search result stats
    try:
        print("trying >1 case")
        # the case where there is something > 1 publications
        search_stat = driver.find_element(By.CLASS_NAME, "publications-count").text
        print(search_stat)
        search_stat_truncated = search_stat.split(" ")[-2]
        print(">1 case done", search_stat_truncated)
    except:
        # handles the case where theres no result
        try:
            print("=0 case")
            search_stat = driver.find_element(By.CSS_SELECTOR,
            'body > app-root > div > app-docsum > div > app-page > div > div.col.l6.m12.s12 > div > div.info-msg > p > b').text
            print(search_stat)
            if "No" in search_stat:
                search_stat_truncated = '0'
                print("=0 case done", search_stat)
        except:
            print("=1 case")
            search_stat_truncated = '0'
    driver.close()
    return search_stat_truncated


def readInputs(path, outputRoot, headless=False):
    """input a csv file with an empty matrix and fill it"""
    df = pd.read_csv(path, encoding='iso-8859-1')
    rowNames = list(df.iloc[:, 0])
    colNames = list(df.columns)[1:]
    outputDict = {}
    outputDict = {list(df.columns)[0]: rowNames}
    outputdf = pd.DataFrame()
    # if the output file already exists, read it in and cancatenate directly
    if os.path.isfile(outputRoot + 'PubtatorOutputs.csv'):
        print("exists")
        existingOutput = pd.read_csv(outputRoot + 'PubtatorOutputs.csv')
        existingCols = list(existingOutput)
        print(existingCols)
        newColNames = [
            colName for colName in colNames if colName not in existingCols]
    else:
        newColNames = colNames

    # iterate through the new columns
    for colName in newColNames:
        outputDict[colName] = []
        for rowName in rowNames:
            meshTerm = '\"' + rowName + '\" AND '
            if '&' not in colName:
                meshTerm += '\"' + colName + '\"'
            else:
                colNameL = colName.split('&')
                meshTerm += '\"' + colNameL[0][:-1] + \
                    '\" AND ' + '\"' + colNameL[1][1:] + '\"'
            print(meshTerm)
            currentNum = pubtatorAutoSearch(
                meshTerm, headless)
            outputDict[colName].append(currentNum)
            print(outputDict[colName])
        outputdf = pd.DataFrame.from_dict(outputDict)

        # if the output file already exists, read it in and cancatenate directly
        if os.path.isfile(outputRoot + 'PubtatorOutputs.csv'):
            existingOutput = pd.read_csv(outputRoot + 'PubtatorOutputs.csv')
            existingCols = list(existingOutput)
            finalOutputdf = pd.concat([existingOutput, outputdf], axis=1)
        else:
            finalOutputdf = outputdf
        finalOutputdf.to_csv(outputRoot + 'PubtatorOutputs.csv', index=False)
        # rest for a bit
        time.sleep(5)
    return outputdf


path = '/Users/fuc/Documents/GitHub/PubtatorSearch/autism_pubtator_source_071924.csv'
outputRoot = '/Users/fuc/Documents/GitHub/PubtatorSearch/outputs/071924_autism_pubtator'
readInputs(path, outputRoot, headless=False)
