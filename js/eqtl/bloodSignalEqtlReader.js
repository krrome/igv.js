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
    igv.BloodSignalEqtlReader = function (config) {

        this.config = config;
        this.url = config.url;
        this.dataset = config.dataset || "merged";
        this.cellType = config.cellType;
        this.indexed = false;
        this.no_subsetting = true;
        this.queryVarId = config.queryVarId;
        this.queryProbeId = config.queryProbeId;
        this.fwdIter = config.fwdIter;
        this.id = config.id;
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

        igv.BloodSignalEqtlReader.prototype.readFeatures = function (chr, bpStart, bpEnd) {
            var self=this

            if (self.url instanceof File) {
                return self.readFeatures_file(chr, bpStart, bpEnd);
            } else {
                return self.readFeatures_rest(chr, bpStart, bpEnd);
            }
        }

        igv.BloodSignalEqtlReader.prototype.readFeatures_rest = function (chr, bpStart, bpEnd) {

            var self=this,
                queryChr = chr.startsWith("chr") ? chr.substr(3) : chr,
                queryStart = Math.floor(bpStart),
                queryEnd = Math.ceil(bpEnd),
                //queryURL = this.url + "all_eqtl?dataset=" + this.dataset + "&cellType=" + this.cellType+"&variantId=" + this.queryVarId + "&probeId=" + this.queryProbeId + "&fwdIter=" + this.fwdIter;
                queryURL = this.url + "all_eqtl?id=" + this.id;

            return new Promise(function (fulfill, reject) {

                igv.xhr.loadJson(queryURL, {
                    withCredentials: self.config.withCredentials
                }).then(function (json) {

                    var variants;

                    if (json) {
                        if (json.detailed_eqtl){
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
                            for (dataset in json.detailed_eqtl){
                                for (queryVarId in json.detailed_eqtl[dataset]){
                                    for (probeId in json.detailed_eqtl[dataset][queryVarId]){
                                        for (fwdIter in json.detailed_eqtl[dataset][queryVarId][probeId]){
                                            ds = json.detailed_eqtl[dataset][queryVarId][probeId][fwdIter]
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
                                                    eqtl.cellType = ds.cellType[i]
                                                    eqtl.dataset = dataset;
                                                    eqtl.queryVarId = queryVarId;
                                                    if (ds.ldR2){
                                                        eqtl.ld = ds.ldR2[i];
                                                    }
                                                    if (eqtl.pValue == 0) {
                                                        eqtl.pValue = 1e-300;
                                                    }
                                                    eqtls.push(eqtl);
                                            }
                                        }
                                    }
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


        igv.BloodSignalEqtlReader.prototype.readFeatures_file = function (chr, start, end) {

            var self = this;

            var options = igv.buildOptions(self.config);    // Add oauth token, if any

            return igv.xhr.loadString(self.config.url, options)
                .then(function (data) {
                    self.header = self.parseHeader(data);
                    return self.parseFeatures(data);   // <= PARSING DONE HERE
                })
                // .catch(function(error) {
                //     console.log(error);
                // })
        };

        igv.BloodSignalEqtlReader.prototype.parseHeader = function (data) {

            var line,
                header,
                dataWrapper;
            var last_line_was_comment = false;
            var delimiter = this.delimiter || "\t";

            dataWrapper = igv.getDataWrapper(data);
            header = {}
            header.skipRows = 0

            while (line = dataWrapper.nextLine()) {
                if (line.startsWith("#")) {
                    if (line.startsWith("#lead")) {
                        header['queryVarId'] = line.split(" ")[1]
                    }
                    else if (line.startsWith("#name")) {
                        if (line.length > 5){
                            header['trackName'] = line.substr(6);
                        }                        
                    }
                    last_line_was_comment = true;
                }
                else {
                    if (last_line_was_comment){
                        // Parse tsv header
                        var columns = {}
                        tokens = line.split(delimiter)
                        for (i in tokens){
                            columns[tokens[i]] = i;
                        }
                        header['columns'] = columns;
                    } else {
                        break;
                    }
                    last_line_was_comment = false;
                }
                header.skipRows++;
            }
            return header;
        };

        igv.BloodSignalEqtlReader.prototype.decode_signal_file = function(tokens) {
            var header = this.header;
            var columns = header.columns;
            get_data_item = function(key){
                if (key in columns){
                    return tokens[columns[key]];
                } else {
                    return null;
                }  
            }
            var eqtl = {};
            chrom = get_data_item("chromosome");
            if ((chrom != null) && (!chrom.startsWith("chr"))) {
                chrom = "chr" + chrom
            }
            eqtl.chr = chrom;
            eqtl.position = get_data_item("start");
            eqtl.start = eqtl.position == null ? eqtl.position : eqtl.position-1;
            eqtl.end = eqtl.position;
            eqtl.snp = get_data_item("snpId");
            eqtl.variantId = get_data_item("variantId");
            eqtl.fwdIter = get_data_item("fwdIter");
            eqtl.geneSymbol = get_data_item("geneSymbol");
            eqtl.beta = get_data_item("beta");
            eqtl.probeId = get_data_item("probeId");
            eqtl.pValue = get_data_item("pValue");
            eqtl.cellType = get_data_item("cellType");
            eqtl.dataset = get_data_item("dataset");
            eqtl.trait = get_data_item("trait");
            queryVarId = get_data_item("queryVarId");
            eqtl.queryVarId = queryVarId == null ? header.queryVarId : queryVarId;
            eqtl.ld = get_data_item("ldR2");
            if (eqtl.pValue == 0) {
                eqtl.pValue = 1e-300;
            }
            return eqtl
        }

        igv.BloodSignalEqtlReader.prototype.parseFeatures = function (data) {

            if (!data) return null;

            var dataWrapper,
                wig,
                feature,
                tokens,
                allFeatures = [],
                line,
                i,
                cnt = 0,
                j,
                decode = this.decode,
                format = this.format,
                delimiter = this.delimiter || "\t";

            dataWrapper = igv.getDataWrapper(data);
            i = -1;

            eqtls = []
            while (line = dataWrapper.nextLine()) {
                i++;
                if (i < this.header.skipRows) continue;

                tokens = line.split(delimiter);
                if (tokens.length != Object.keys(this.header.columns).length) {
                    continue;
                }

                //token assignment here:
                eqtls.push(this.decode_signal_file(tokens));
            }

            return eqtls;
        };


    return igv;
})
(igv || {});
