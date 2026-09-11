#target illustrator
(function(){
  app.userInteractionLevel=UserInteractionLevel.DONTDISPLAYALERTS;
  var b="D:/桌面/sepsis/文章/返修投稿新主图/";
  var d=app.open(new File("D:/桌面/sepsis/review_revision/09_figures/S8_viability.svg"));
  var t=new ExportOptionsTIFF();t.imageColorSpace=ImageColorSpace.RGB;
  t.resolution=300;t.antiAliasing=AntiAliasingMethod.ARTOPTIMIZED;
  t.lZWCompression=true;t.byteOrder=TIFFByteOrder.IBMPC;t.embedICCProfile=true;
  t.saveMultipleArtboards=true;t.artboardRange="1";
  var tf=new File(b+"S8_fig.tif");d.exportFile(tf,ExportType.TIFF,t);
  var alt=new File(b+"S8_fig-01.tif");if(alt.exists){if(tf.exists)tf.remove();alt.rename("S8_fig.tif");}
  var po=new PDFSaveOptions();po.compatibility=PDFCompatibility.ACROBAT7;
  po.preserveEditability=true;po.optimization=true;
  d.saveAs(new File(b+"S8_fig.pdf"),po);
  var ao=new IllustratorSaveOptions();ao.pdfCompatible=true;ao.compressed=true;
  d.saveAs(new File(b+"S8_fig.ai"),ao);d.close(SaveOptions.DONOTSAVECHANGES);
  var m=new File(b+"S8_fig_python_master.pdf");if(m.exists)m.remove();
  app.userInteractionLevel=UserInteractionLevel.DISPLAYALERTS;
})();
