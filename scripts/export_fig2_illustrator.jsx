#target illustrator

(function () {
    var masterFile = new File("D:/桌面/sepsis/文章/返修投稿新主图/Fig2_python_master.pdf");
    var aiFile = new File("D:/桌面/sepsis/文章/返修投稿新主图/Fig2.ai");
    var pdfFile = new File("D:/桌面/sepsis/文章/返修投稿新主图/Fig2.pdf");
    var tifFile = new File("D:/桌面/sepsis/文章/返修投稿新主图/Fig2.tif");

    if (!masterFile.exists) {
        throw new Error("Fig2 Python vector master not found: " + masterFile.fsName);
    }

    app.userInteractionLevel = UserInteractionLevel.DONTDISPLAYALERTS;
    var doc = app.open(masterFile);

    // Keep every label live and editable.  Apply the journal-safe Arial face
    // text frame by text frame without converting type to outlines.
    var arial = null;
    var arialBold = null;
    try {
        arial = app.textFonts.getByName("ArialMT");
    } catch (fontError) {
        try {
            arial = app.textFonts.getByName("Arial");
        } catch (fallbackError) {
            arial = null;
        }
    }
    try {
        arialBold = app.textFonts.getByName("Arial-BoldMT");
    } catch (boldFontError) {
        arialBold = arial;
    }
    if (arial !== null) {
        for (var i = 0; i < doc.textFrames.length; i++) {
            try {
                var frame = doc.textFrames[i];
                var content = frame.contents.replace(/^\s+|\s+$/g, "");
                var attr = frame.textRange.characterAttributes;
                var oldName = "";
                try { oldName = attr.textFont.name; } catch (nameError) {}
                var keepBold = /Bold/i.test(oldName) || /^[A-F]$/.test(content);
                attr.textFont = keepBold ? arialBold : arial;
                if (/^[A-F]$/.test(content)) { attr.size = 12; }
            } catch (frameError) {}
        }
    }

    // RGB is appropriate for PLOS online figures and avoids unintended CMYK
    // conversion of the restrained colour palette.
    try {
        doc.documentColorSpace = DocumentColorSpace.RGB;
    } catch (colorSpaceError) {}

    var tifOptions = new ExportOptionsTIFF();
    tifOptions.imageColorSpace = ImageColorSpace.RGB;
    tifOptions.resolution = 300;
    tifOptions.antiAliasing = AntiAliasingMethod.ARTOPTIMIZED;
    tifOptions.lZWCompression = true;
    tifOptions.byteOrder = TIFFByteOrder.IBMPC;
    tifOptions.embedICCProfile = true;
    tifOptions.saveMultipleArtboards = true;
    tifOptions.artboardRange = "1";
    doc.exportFile(tifFile, ExportType.TIFF, tifOptions);
    // Illustrator appends the artboard suffix when artboard clipping is used.
    // Normalize it back to the manuscript's stable Fig2.tif filename.
    var artboardTif = new File("D:/桌面/sepsis/文章/返修投稿新主图/Fig2-01.tif");
    if (artboardTif.exists) {
        if (tifFile.exists) {
            tifFile.remove();
        }
        artboardTif.rename("Fig2.tif");
    }

    var aiOptions = new IllustratorSaveOptions();
    aiOptions.pdfCompatible = true;
    aiOptions.embedICCProfile = true;
    aiOptions.compressed = true;
    doc.saveAs(aiFile, aiOptions);

    // Save PDF last. Illustrator can otherwise retain an existing PDF when a
    // subsequent AI save changes the active document target.
    var pdfOptions = new PDFSaveOptions();
    pdfOptions.compatibility = PDFCompatibility.ACROBAT7;
    pdfOptions.preserveEditability = true;
    pdfOptions.generateThumbnails = true;
    pdfOptions.optimization = true;
    pdfOptions.viewAfterSaving = false;
    if (pdfFile.exists) { pdfFile.remove(); }
    doc.saveAs(pdfFile, pdfOptions);

    doc.close(SaveOptions.DONOTSAVECHANGES);
    // The submission folder keeps only the three requested deliverables.
    var transientFiles = [
        new File("D:/桌面/sepsis/文章/返修投稿新主图/Fig2.png"),
        new File("D:/桌面/sepsis/文章/返修投稿新主图/Fig2.svg"),
        new File("D:/桌面/sepsis/文章/返修投稿新主图/Fig2_preview.png"),
        new File("D:/桌面/sepsis/文章/返修投稿新主图/Fig2_illustrator_preview.png"),
        new File("D:/桌面/sepsis/文章/返修投稿新主图/Fig2_python_master.pdf"),
        new File("D:/桌面/sepsis/review_revision/scripts/07_class_imbalance_summary.R")
    ];
    for (var j = 0; j < transientFiles.length; j++) {
        if (transientFiles[j].exists) {
            transientFiles[j].remove();
        }
    }
    app.userInteractionLevel = UserInteractionLevel.DISPLAYALERTS;
})();
