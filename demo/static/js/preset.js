document.getElementById('WGSBtn').addEventListener('click', function() {initialPreset("WGS",fq_type);});


async function initialPreset(presetType, fq_type) {
    try {
        const response = await fetch('/batch_code/get_preset', {
            method: 'POST', headers: {'Content-Type': 'application/json'},
            body: JSON.stringify({presetType: presetType,  mode: fq_type})});
        const data = await response.json();  // 这里返回soft_step列表
       

 
        data.options.forEach(element => {addStep(element.type)});
        initOptionsUI(data)
     	disableAllButtons()
    } catch (error) { console.error('Error initialPreset:', error);}
}

function disableAllButtons() {
    var presetIds = ["WGSBtn", "WGBSBtn", "CSBtn", "CAPSBtn"];
    var buttonIds = [];
    Object.keys(available_steps).forEach(function(category) {
        available_steps[category].forEach(function(button) {
            buttonIds.push(`btn-${category}-${button}`); });
    });
    var allButtonIds = buttonIds.concat(presetIds);
    allButtonIds.forEach(function(id) {
        if (document.getElementById(id)) {document.getElementById(id).disabled = true;}});
}


