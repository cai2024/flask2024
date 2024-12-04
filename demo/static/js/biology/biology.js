$(document).ready(function() {
    initialCustom();
    document.getElementById('confirmBtn').addEventListener('click', confirmSteps);
    document.getElementById('paramBtn').addEventListener('click', paramSteps);
    document.getElementById('submitBtn').addEventListener('click', submitSteps);
});


async function initialCustom() {
    const presetKey = `preset_${form["pipeline"]}_${form["mode"]}`;   
    if (form["param_choice"] === "default") {initialPreset(presets_dict[presetKey]); return;}
    const responseData = await postData('/biology/get_preset', {presetKey: presets_dict[presetKey]});
    steps = responseData.steps;
    responseData.soft_list.forEach((stepObj, index) => {
        const step = Object.keys(stepObj)[0];
        let html = `<div id="input-${step}" class="form-group">
                        <label><strong>${step}:</strong></label>
                        <select id="select-${step}" class="form-control mt-2">`;
        stepObj[step].forEach(option => {html += `<option value="${option}">${option}</option>`;});
        html += `</div>`;
        document.getElementById('selected-steps').insertAdjacentHTML('beforeend', html);
    });
}


async function confirmSteps() {
    document.getElementById('selected-steps').querySelectorAll('.form-group').forEach((step,i) => {
        const select = step.querySelector('select');
        steps[i].type = select.options[select.selectedIndex].value });  
    const responseData = await postData('/batch_code/init_div_steps', {steps});
    if (responseData) {
        steps = responseData
        initOptionsUI(responseData);
        confirmBtn.disabled = true;
        document.querySelectorAll('#selected-steps select').forEach(select => {select.disabled = true;})
        document.getElementById('paramBtn').disabled = false;
    } 
}



async function paramSteps() {
    try {
        collectStepInputs(steps);
        const responseData = await postData('/biology/submit_config_script', {steps});
        document.getElementById('dynamicCode').textContent = responseData.my_code;
        hljs.highlightElement(document.getElementById('dynamicCode'));
        document.getElementById('submitBtn').disabled = false;
        document.querySelectorAll('#inputArea input, #inputArea select, #indexArea select').forEach(el => el.disabled = true);
    } catch (error) {console.error('Error submitting steps:', error);}
}



async function submitSteps() {
    try {
        const response = await fetch('/common/submit_smk_job', {
            method: 'POST', headers: {'Content-Type': 'application/json'},
            body: JSON.stringify({ task_id: task_id })
        }); if (response.ok) {window.location.href = '/query/' + task_id;} 
    } catch (error) {console.error('Error submitting job:', error);}
}
