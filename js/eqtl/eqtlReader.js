/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2015 UC San Diego
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

/**
 * Created by jrobinso on 10/8/15.
 */

var igv = (function (igv) {


    /**
     * @param url - url to the webservice
     * @constructor
     */
    igv.BloodEqtlReader = function (config) {

        this.config = config;
        this.url = config.url;
        this.dataset = config.dataset;
        this.indexed = true;
    };

    //{
    //    "release": "v6",
    //    "singleTissueEqtl": [
    //    {
    //        "beta": -0.171944779728988,
    //        "chromosome": "3",
    //        "gencodeId": "ENSG00000168827.10",
    //        "geneSymbol": "GFM1",
    //        "pValue": 1.22963421134407e-09,
    //        "snpId": "rs3765025",
    //        "start": 158310846,
    //        "tissueName": "Thyroid"
    //    },
    //
    // http://vgtxportaltest.broadinstitute.org:9000/v6/singleTissueEqtlByLocation?tissueName=Thyroid&chromosome=3&start=158310650&end=158311650

        igv.BloodEqtlReader.prototype.readFeatures = function (chr, bpStart, bpEnd) {

            var self=this,
                queryChr = chr.startsWith("chr") ? chr.substr(3) : chr,
                queryStart = Math.floor(bpStart),
                queryEnd = Math.ceil(bpEnd),
                //TODO: enable querying the high-detail variant only if the region is within 1.5MB of the specified variant.
                queryURL = this.url + "eqtl/" + this.dataset + "/" + queryChr + "/" + queryStart + "/" + queryEnd;

            return new Promise(function (fulfill, reject) {

                igv.xhr.loadJson(queryURL, {
                    withCredentials: self.config.withCredentials
                }).then(function (json) {

                    var variants;

                    if (json) {
                        if (json.lead_eqtl){
                            //variants = json.variants;
                            //variants.sort(function (a, b) {
                            //    return a.POS - b.POS;
                            //});
                            //source.cache = new FeatureCache(chr, queryStart, queryEnd, variants);

                            /*
                                _eqtl_data_ = {
                                   "beta": [0.1,0.2,0.3],
                                   "chromosome": ["12", "12", "12"],
                                   "geneSymbol": ["APAF1", "IKBIP", "XYZ"],
                                   "probeId": ["23423", "24524235", "2323423"],
                                   "pValue": [1.37476e-10, 1.37476e-10, 1.37476e-10], 
                                   "snpId": ["rs111604160", "rs111604160", "rs111604160"], 
                                   "start": [98980788, 98980788, 98980788],
                                   "variantId": ["12:98980788:A:G", "12:98980788:A:G", "12:98980788:A:G"],
                                   "fwdIter": [0,1,0],
                                   "ldR2": None / [.8, .9, .2],
                                }
                            */

                            var eqtls = [];

                            //TODO: keep the structure with cell type and also return the high-detail eQTL results!
                            for (dataset in json.lead_eqtl){
                                ds = json.lead_eqtl[dataset]
                                ds_eqtls = []
                                for (i =0; i < ds.chromosome.length; i++){
                                    var eqtl = {};
                                    eqtl.chr = "chr" + ds.chromosome[i];
                                    eqtl.position = ds.start[i];
                                    eqtl.start = ds.start[i]-1;
                                    eqtl.end = ds.start[i];
                                    eqtl.snp = ds.snpId[i];
                                    eqtl.variantId = ds.variantId[i];
                                    eqtl.fwdIter = ds.fwdIter[i];
                                    eqtl.geneSymbol = ds.geneSymbol[i];
                                    eqtl.beta = ds.beta[i];
                                    eqtl.probeId = ds.probeId[i];
                                    eqtl.pValue = ds.pValue[i];
                                    eqtl.dataset = dataset
                                    eqtls.push(eqtl);
                                }
                            }

                            fulfill(eqtls);

                        }
                        else {
                            fulfill(null);
                        }
                    }
                    else {
                        fulfill(null);
                    }

                }).catch(reject);

            });
        }


    return igv;
})
(igv || {});
