/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2014 Broad Institute
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

var igv = (function (igv) {


    igv.BloodEqtlTrack = function (config) {

        /*
        eqtl.chr = "chr" + json.lead_eqtl.chromosome[i];
        eqtl.position = json.lead_eqtl.start[i];
        eqtl.start = json.lead_eqtl.start[i]-1;
        eqtl.end = json.lead_eqtl.start[i];
        eqtl.snp = json.lead_eqtl.snpId[i];
        eqtl.variantId = json.lead_eqtl.variantId[i];
        eqtl.fwdIter = json.lead_eqtl.fwdIter[i];
        eqtl.geneSymbol = json.lead_eqtl.geneSymbol[i];
        eqtl.beta = json.lead_eqtl.beta[i];
        eqtl.probeId = json.lead_eqtl.probeId[i];
        eqtl.pValue = json.lead_eqtl.pValue[i];
        eqtls.push(eqtl);
        */


        var url = config.url,
            label = config.name;

        this.config = config;
        this.url = url;
        this.name = label;
        this.dataset = config.dataset || "merged";
        this.pValueField = config.pValueField || "pValue";
        this.geneField = config.geneField || "geneSymbol";
        this.snpField = config.snpField || "snp";

        this.autoscale = (config.autoScale === undefined ? true : config.autoScale);
        this.percentile = (config.percentile === undefined ? 98 : config.percentile);

        this.minLogP = config.minLogP || 3.5;
        this.maxLogP = config.maxLogP || 25;
        this.background = config.background;    // No default
        this.divider = config.divider || "rgb(225,225,225)";
        this.dotSize = config.dotSize || 2;
        this.height = config.height || 100;
        this.autoHeight = false;
        this.disableButtons = config.disableButtons;
        this.visibilityWindow = config.visibilityWindow;

        this.featureSource = new igv.FeatureSource(config);

    };

    igv.BloodEqtlTrack.prototype.paintAxis = function (ctx, pixelWidth, pixelHeight) {

        var track = this,
            yScale = (track.maxLogP - track.minLogP) / pixelHeight;

        var font = {
            'font': 'normal 10px Arial',
            'textAlign': 'right',
            'strokeStyle': "black"
        };

        igv.graphics.fillRect(ctx, 0, 0, pixelWidth, pixelHeight, {'fillStyle': "rgb(255, 255, 255)"});
        
        // Determine a tick spacing such that there is at least 10 pixels between ticks
        
        var n = Math.ceil((this.maxLogP - this.minLogP) * 10 / pixelHeight);
        

        for (var p = 4; p <= track.maxLogP; p += n) {

            var x1,
                x2,
                y1,
                y2,
                ref;

            // TODO: Dashes may not actually line up with correct scale. Ask Jim about this

            ref = 0.85 * pixelWidth;
            x1 = ref - 5;
            x2 = ref;

            y1 = y2 = pixelHeight - Math.round((p - track.minLogP) / yScale);

            igv.graphics.strokeLine(ctx, x1, y1, x2, y2, font); // Offset dashes up by 2 pixel

            if(y1 > 8) {
                igv.graphics.fillText(ctx, p, x1 - 1, y1 + 2, font);
            } // Offset numbers down by 2 pixels;
        }

        font['textAlign'] = 'center';

        igv.graphics.fillText(ctx, "-log10(pvalue)", pixelWidth / 4, pixelHeight / 2, font, {rotate: {angle: -90}});

    };

    igv.BloodEqtlTrack.prototype.getFeatures = function (chr, bpStart, bpEnd) {
        return this.featureSource.getFeatures(chr, bpStart, bpEnd);
    };

    igv.BloodEqtlTrack.prototype.draw = function (options) {

        var self = this,
            featureList = options.features,
            ctx = options.context,
            bpPerPixel = options.bpPerPixel,
            bpStart = options.bpStart,
            pixelWidth = options.pixelWidth,
            pixelHeight = options.pixelHeight,
            bpEnd = bpStart + pixelWidth * bpPerPixel + 1,
            yScale = (self.maxLogP - self.minLogP) / pixelHeight,
            selection = options.genomicState.selection;

        /*
        // make sure that only eQTLs are plotted that belong to the given dataset.
        var featureList_in_dataset = []
        for (el in featureList) {
            if (featureList[el].dataset == this.config.dataset) featureList_in_dataset.push(featureList[el]);
        }
        featureList = featureList_in_dataset
        */

        // Background
        if (this.background) igv.graphics.fillRect(ctx, 0, 0, pixelWidth, pixelHeight, {'fillStyle': this.background});
        igv.graphics.strokeLine(ctx, 0, pixelHeight - 1, pixelWidth, pixelHeight - 1, {'strokeStyle': this.divider});

        if (ctx) {

            var len = featureList.length;
            
            ctx.save();
            
            self.maxLogP = autoscale(featureList, bpStart, bpEnd);
            yScale = (self.maxLogP - self.minLogP) / pixelHeight;
            
            // Draw in two passes, with "selected" eqtls drawn last
            drawEqtls(false);
            drawEqtls(true);
            
            ctx.restore();

        }

        function autoscale(featureList, start, end) {

            var values = featureList
                .filter(function (eqtl) {
                    return eqtl.position >= start && eqtl.position <= end
                })
                .map(function (eqtl) {
                    return -Math.log(eqtl[self.pValueField]) / Math.LN10
                })

            return igv.Math.percentile(values, self.percentile);
            
        }

        function drawEqtls(drawSelected) {

            var radius = drawSelected ? 2 * self.dotSize : self.dotSize,
                eqtl,
                i,
                px,
                py,
                color,
                isSelected,
                snp,
                geneName,
                capped;

            for (i = 0; i < len; i++) {

                eqtl = featureList[i];
                px = (Math.round(eqtl.position - bpStart + 0.5)) / bpPerPixel;
                if (px < 0) continue;
                else if (px > pixelWidth) break;


                snp = eqtl.snp.toUpperCase();
                geneName = eqtl[self.geneField].toUpperCase();

                isSelected = selection &&
                    (selection.snp === snp || selection.gene === geneName);

                if (!drawSelected || isSelected) {

                    // Add eqtl's gene to the selection if this is the selected snp.
                    // TODO -- this should not be done here in the rendering code.
                    if (selection && selection.snp === snp) {
                        selection.addGene(geneName);
                    }

                    var mLogP = -Math.log(eqtl[self.pValueField]) / Math.LN10;
                    if (mLogP >= self.minLogP) {

                        if(mLogP > self.maxLogP) {
                            mLogP = self.maxLogP;
                            capped = true;
                        } else {
                            capped = false;

                        }

                        py = Math.max(0 + radius, pixelHeight - Math.round((mLogP - self.minLogP) / yScale));
                        eqtl.px = px;
                        eqtl.py = py;

                        if (drawSelected && selection) {
                            color = selection.colorForGene(geneName);
                            igv.graphics.setProperties(ctx, {fillStyle: color, strokeStyle: "black"});
                        } else {
                            color = capped ? "rgb(10, 10, 10)" : "rgb(60, 60, 60)";
                            igv.graphics.setProperties(ctx, {fillStyle: color, strokeStyle: color});
                        }

                        igv.graphics.fillCircle(ctx, px, py, radius);
                        igv.graphics.strokeCircle(ctx, px, py, radius);
                    }
                }
            }
        }

    };

    igv.BloodEqtlTrack.prototype.menuItemList = function () {
        return [];
    };

    /**
     * Return "popup data" for feature @ genomic location.  Data is an array of key-value pairs
     */
    igv.BloodEqtlTrack.prototype.popupData = function (config) {
    
        var genomicLocation = config.genomicLocation,
            xOffset = config.x,
            yOffset = config.y,
            referenceFrame = config.viewport.genomicState.referenceFrame;

        // We use the featureCache property rather than method to avoid async load.  If the
        // feature is not already loaded this won't work,  but the user wouldn't be mousing over it either.
        if (config.viewport.tile.features) {

            var chr = referenceFrame.chrName,
                tolerance = 2 * this.dotSize * referenceFrame.bpPerPixel,
                featureList = config.viewport.tile.features,
                dotSize = this.dotSize,
                tissue = this.name;

            if (featureList && featureList.length > 0) {


                var popupData = [];
                featureList.forEach(function (feature) {
                    if (feature.end >= genomicLocation - tolerance &&
                        feature.start <= genomicLocation + tolerance &&
                        feature.py - yOffset < 2 * dotSize) {

                        if (popupData.length > 0) {
                            popupData.push("<hr>");
                        }

                        detailsTrackName = feature.dataset + " eQTL " + feature.geneSymbol + " (" + feature.probeId + ") with var " + feature.variantId + " ("+feature.fwdIter+")";

                        popupData.push(
                            {name: "cell type", value: feature.cellType},
                            {name: "gene hgnc", value: feature.geneSymbol},
                            {name: "variant id", value: feature.variantId},
                            {name: "rs-id", value: feature.snp},
                            {name: "p value", value: feature.pValue},
                            {name: "effect size", value: feature.beta},
                            {name: "probe id", value: feature.probeId},
                            {name: "detailed", value: "Load details", href: "create_track_by_id('"+feature.id+"', hgnc='"+feature.geneSymbol+"')"}
                            );

                    }
                });
                return popupData;
            }
        }
    }
    

    return igv;

})(igv || {});
