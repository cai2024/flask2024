document.getElementById('submitBtn').addEventListener('click', submitSteps);
document.getElementById('unlockBtn').addEventListener('click', function() {initialPreset("WGS");});




// 上传fq_type, raw_myOptions[0], steps到session




async function submitSteps() {
    try {
        const response = await fetch('/submit_form_steps', {
            method: 'POST', headers: {'Content-Type': 'application/json'},
            body: JSON.stringify({fq_type: raw_mode, raw_myOptions: raw_myOptions[0], steps: steps})
        });
        if (response.ok) {window.location.href = '/submit_form_code/' + task_id;}}
    catch (error) {console.error('Error submitting steps:', error);}
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
