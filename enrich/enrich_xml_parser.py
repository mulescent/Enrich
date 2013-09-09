#!/usr/bin/env python
'''
enrich_xml_parser: this module implements the SAX XML parser to parse enrich xml configuration files
'''

import sys, optparse, xml.sax #imports of standard modules

def main(infile):
    
    #results are written to by the enrichXMLHandler
    results = {}
    handler = enrichXMLHandler(results)
    
    #create a SAX XML parser instance
    parser = xml.sax.make_parser()
    parser.setContentHandler(handler)
    
    #parse the input file
    try:
        parser.parse(infile)
        
    except xml.sax.SAXParseException:
        sys.exit('Error: improperly formatted XML file')
    
    #return results 
    if __name__ == '__main__':
        print results
    
    else:
        return(results)
    
class enrichXMLHandler(xml.sax.ContentHandler):
    '''enrichXMLHandler: this is a SAX ContentHandler class for enrich configuration files'''
    
    def __init__(self, results):
        self.results = results
        self.sectName = '' #stores the current section
        self.inSect = '' #stores the current element's section
        self.elemName = '' #stores the current element's name
        self.isRoot = int() #stores the current element's root_element attribute
        self.isSection = int() #stores the current element's section attribute
        
    def startElement(self, name, attrs):

        if attrs.get('root_element') == 'TRUE':
            self.isRoot = 1
            
        elif attrs.get('section') == 'TRUE':
            self.isSection = 1
            self.isRoot = 0
            self.results[str(name)] = {} #must use str() function to convert name to plain text from unicode, xml.sax handlers always return unicode
            self.sectName = str(name) 
            self.elemName = str(name)
        
        else:
            self.isRoot = 0
            self.isSection = 0
            self.elemName = str(name)
            
        self.inSect = self.sectName
        
    def characters(self, ch):
    
        if self.isSection == 0 and self.isRoot != 1: #check to make sure this character is associated with neither a root or section element
            
            if self.elemName not in self.results[self.inSect]: #do not overwrite an existing value
            
                self.results[self.inSect][self.elemName] = str(ch).rstrip()             
            
if __name__ == '__main__':
    parser = optparse.OptionParser() #initialize parser to grab command line parameters
    parser.add_option('--path', action = 'store', type = 'string', dest = 'path', help = 'path to project directory')
    parser.add_option('--infile', action = 'store', type = 'string', dest = 'infile', help = 'name of input file')
    option, args = parser.parse_args()
     
    main(option.path + option.infile) 