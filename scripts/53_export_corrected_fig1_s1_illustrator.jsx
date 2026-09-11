#target illustrator

(function () {
    app.userInteractionLevel = UserInteractionLevel.DONTDISPLAYALERTS;

    var base = "C:/Users/ZCL/sepsis_final_figures/";
    var regular = null;
    var bold = null;
    try { regular = app.textFonts.getByName("ArialMT"); } catch (e1) {}
    try { bold = app.textFonts.getByName("Arial-BoldMT"); } catch (e2) {}

    function normalizeFonts(doc) {
        for (var i = 0; i < doc.textFrames.length; i++) {
            var tf = doc.textFrames[i];
            var oldName = "";
            try { oldName = tf.textRange.characterAttributes.textFont.name.toLowerCase(); } catch (e3) {}
            var isBold = oldName.indexOf("bold") >= 0 || oldName.indexOf("black") >= 0;
            if (isBold && bold !== null) {
                tf.textRange.characterAttributes.textFont = bold;
            } else if (regular !== null) {
                tf.textRange.characterAttributes.textFont = regular;
            }
        }
    }

    function exportFigure(sourceName, outputStem) {
        var doc = app.open(new File(base + sourceName));
        normalizeFonts(doc);

        var tiffOptions = new ExportOptionsTIFF();
        tiffOptions.imageColorSpace = ImageColorSpace.RGB;
        tiffOptions.resolution = 300;
        tiffOptions.antiAliasing = AntiAliasingMethod.ARTOPTIMIZED;
        tiffOptions.lZWCompression = true;
        tiffOptions.byteOrder = TIFFByteOrder.IBMPC;
        tiffOptions.embedICCProfile = true;
        tiffOptions.saveMultipleArtboards = true;
        tiffOptions.artboardRange = "1";
        var tif = new File(base + outputStem + ".tif");
        doc.exportFile(tif, ExportType.TIFF, tiffOptions);
        var artboardTif = new File(base + outputStem + "-01.tif");
        if (artboardTif.exists) {
            if (tif.exists) tif.remove();
            artboardTif.rename(outputStem + ".tif");
        }

        var pdfOptions = new PDFSaveOptions();
        pdfOptions.compatibility = PDFCompatibility.ACROBAT7;
        pdfOptions.preserveEditability = true;
        pdfOptions.optimization = true;
        pdfOptions.viewAfterSaving = false;
        doc.saveAs(new File(base + outputStem + ".pdf"), pdfOptions);

        var aiOptions = new IllustratorSaveOptions();
        aiOptions.pdfCompatible = true;
        aiOptions.embedICCProfile = true;
        aiOptions.compressed = true;
        doc.saveAs(new File(base + outputStem + ".ai"), aiOptions);
        doc.close(SaveOptions.DONOTSAVECHANGES);
    }

    exportFigure("Fig1_corrected_candidate.pdf", "Fig1_corrected_final");
    exportFigure("S1_fig_R_master.pdf", "S1_fig_corrected_final");
    app.userInteractionLevel = UserInteractionLevel.DISPLAYALERTS;
})();
