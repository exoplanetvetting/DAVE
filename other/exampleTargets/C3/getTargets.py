# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 09:45:08 2015

@author: fmullall
"""


import matplotlib.pyplot as mp
import numpy as np

import queryMast
import nca

class QueryK2ByTgtId(queryMast.QueryMast):
    def __init__(self):
        """Query the MAST for K2 data by Investigation ID

        All the code that does the work is stored in QueryMast and
        AbstractDataQuery. This class just needs to specifiy the
        name of the KIC at MAST.
        """

        self.cachePrefix = "k2invest"
        queryMast.QueryMast.__init__(self)
        self.catalogueName = "k2/data_search"
        self.cachePrefix = "k2invest"


    def constructUrl(self, investigationId):

        if self.catalogueName is None:
            raise NotImplementedError("Daughter class must set self.catalogueName!")

        #Cast as string if necessary
        if isinstance(investigationId, int):
            investigationId = "%i" %(investigationId)

        outputCols = "ktc_k2_id,sci_campaign,ktc_investigation_id,sci_ra,sci_dec,kp"

        # build mast query
        url  = 'http://archive.stsci.edu/'
        url += '%s/search.php?' %(self.catalogueName)
        url += 'action=Search'
        url += "&ktc_investigation_id=*%s*" %(investigationId)
        url += '&coordformat=dec'
        url += '&outputformat=CSV'
        url += '&selectedColumnsCsv=%s' %(outputCols)
        url += '&verb=0'

        return url


    def queryProposal(self, investigationId):
        url = self.constructUrl(investigationId)
        cacheDir = "cache"

        text = self.query(0,0,0, url=url, cacheDir=cacheDir)
        return text


def main():
    props = []
    props.append( (3111, 'WD_Kilic'))
    props.append( (3116, 'WD_Redfield'))
    props.append( (3096, 'EXO_Heller'))
    props.append( (3095, 'EXO_vanGrootel'))
    props.append( (3086, 'EB_Southworth'))
    props.append( (3067, 'EB_Peters'))
    props.append( (3049, 'EB_Prsa'))
    props.append( (3005, 'EB_Shporer'))
#
    combined = []
    qm = QueryK2ByTgtId()
    for p in props:
        print "Getting %s" %(p[1])
        text = qm.queryProposal(p[0])
        data = qm.parseResults(text)
        data = appendCol(data[1:, :], 'Proposal')
        data[:, 'Proposal'] = p[1]
        combined.append(data)

    cat = np.concatenate(combined)
    idx=  np.unique(cat[:,0], return_index=True)[1]
    cat = cat[idx]
    cat = nca.Nca(cat)
    cat.setLookup(1, data.lookup[1])

    return cat


import tools
def printCat(cat):
    text = []
    for row in cat:
        line = " ".join(row)
        text.append(line)


    headerStr="Example K2C3 stars for algortihm testing"
    colNames = cat.lookup[1]
    hdr = tools.createHeader(headerStr, columnNames=colNames)

    fp = open("K2C3cat.txt", 'w')
    fp.write("\n".join(hdr))
    fp.write("\n")
    fp.write("\n".join(tools.respace(text)))
    fp.close()


def appendCol(data, name):
    nr, nc = data.shape
    nc +=1

    newCol = np.atleast_2d(data[:,0].asarray()).transpose()
    newData = np.hstack( (data.asarray(), newCol))
    lookup = data.lookup
    lookup[1].append(name)

    return nca.Nca(newData, lookup)
