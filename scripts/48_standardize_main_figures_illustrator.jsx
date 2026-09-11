#target illustrator
(function(){
  app.userInteractionLevel=UserInteractionLevel.DONTDISPLAYALERTS;
  var base="D:/桌面/sepsis/文章/返修投稿新主图/";
  var regular=null,bold=null;
  try{regular=app.textFonts.getByName("ArialMT");}catch(e1){try{regular=app.textFonts.getByName("Arial");}catch(e2){}}
  try{bold=app.textFonts.getByName("Arial-BoldMT");}catch(e3){bold=regular;}
  function isSinglePanelLetter(s){return /^[A-H]$/.test(s.replace(/^\s+|\s+$/g,""));}
  var regularPanelSites={
    4:{"A":[1.5,596.3],"B":[265.0,596.3],"C":[1.5,394.3],"D":[62.0,190.3]},
    5:{"A":[1.5,600.3],"B":[265.0,600.3],"C":[1.5,450.3],"D":[265.0,450.3],
       "E":[1.5,300.3],"F":[265.0,300.3],"G":[1.5,150.3],"H":[265.0,150.3]}
  };
  function isKnownRegularPanelSite(n,s,tf){
    var table=regularPanelSites[n], key=s.replace(/^\s+|\s+$/g,"");
    if(!table||!table[key])return false;
    var p=table[key]; return Math.abs(tf.left-p[0])<2.0&&Math.abs(tf.top-p[1])<3.0;
  }
  for(var n=1;n<=8;n++){
    var stem="Fig"+n, pdf=new File(base+stem+".pdf"), doc=app.open(pdf);
    for(var i=0;i<doc.textFrames.length;i++){
      var tf=doc.textFrames[i], content="", attr=null;
      try{content=tf.contents; attr=tf.textRange.characterAttributes;}catch(frameError){continue;}
      var oldName=""; try{oldName=attr.textFont.name;}catch(e4){}
      // A single A-H is a panel label only when it was already bold in the
      // assembled PDF. This prevents experimental group ticks such as C from
      // being mistaken for panel letters and enlarged to 12 pt.
      var panelLabel=isSinglePanelLetter(content)&&
        (/Bold/i.test(oldName)||isKnownRegularPanelSite(n,content,tf));
      try{if(regular!==null) attr.textFont=(/Bold/i.test(oldName)||panelLabel)?bold:regular;}catch(fontError){}
      if(panelLabel){
        try{attr.size=12;}catch(sizeError){}
      }else{
        // PLOS ONE requires text within figures to be 8-12 pt at the final
        // submitted dimensions.  Apply the limit after panel assembly, since
        // half-width placement scales otherwise adequate source-panel type.
        try{
          if(attr.size<8){attr.size=8;}
          else if(attr.size>12){attr.size=12;}
        }catch(finalSizeError){}
      }
    }
    var po=new PDFSaveOptions(); po.compatibility=PDFCompatibility.ACROBAT7; po.preserveEditability=true; po.optimization=true; po.viewAfterSaving=false;
    doc.saveAs(new File(base+stem+".pdf"),po);
    var ao=new IllustratorSaveOptions(); ao.pdfCompatible=true; ao.embedICCProfile=true; ao.compressed=true;
    doc.saveAs(new File(base+stem+".ai"),ao);
    var to=new ExportOptionsTIFF(); to.imageColorSpace=ImageColorSpace.RGB; to.resolution=300; to.antiAliasing=AntiAliasingMethod.ARTOPTIMIZED; to.lZWCompression=true; to.byteOrder=TIFFByteOrder.IBMPC; to.embedICCProfile=true; to.saveMultipleArtboards=true; to.artboardRange="1";
    var tfout=new File(base+stem+".tif"); doc.exportFile(tfout,ExportType.TIFF,to);
    var alt=new File(base+stem+"-01.tif"); if(alt.exists){if(tfout.exists)tfout.remove();alt.rename(stem+".tif");}
    doc.close(SaveOptions.DONOTSAVECHANGES);
  }
  app.userInteractionLevel=UserInteractionLevel.DISPLAYALERTS;
})();
