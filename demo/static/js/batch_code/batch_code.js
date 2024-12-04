// // 审阅时间：20241110
$(document).ready(function() {
    initializeStepButtons('.step-button', 'input[name="fqType"]:checked', 'input[name="fileType"]:checked');
    document.getElementById('confirmBtn').addEventListener('click', confirmSteps);
    document.getElementById('submitBtn').addEventListener('click', submitSteps);
    document.getElementById('lookBtn').addEventListener('click', lookSteps);
    document.getElementById('presetConfirm').addEventListener('click', function() {
        initialPreset(document.getElementById('presetInput').value);
        this.disabled = true;});

    document.getElementById('WGSBtn').addEventListener('click', function() {initialPreset(presets_dict["preset_WGS_pair"])});


});




async function lookSteps() {
    try {
        collectStepInputs(steps);
        const responseData = await postData('/batch_code/display_div_steps', {steps});
        document.getElementById('outputArea').innerText = responseData.base64_encoded_data;
        outputArea.appendChild(document.createElement('br'));
        document.getElementById('outputArea').appendChild(Object.assign(document.createElement('img'),{ src: responseData.image, alt: '生成的图像' }));
        } catch (error) {console.error('Error:', error);
        document.getElementById('outputArea').innerText = '获取数据失败，请重试。';}
    }

async function confirmSteps() {
    for (let i = 0; i < steps.length; i++) { steps[i].id = i; }
    const responseData = await postData('/batch_code/init_div_steps', {steps});
    if (responseData) {
        steps = responseData;
        initOptionsUI(responseData);
        confirmBtn.disabled = true;
        submitBtn.disabled = false;
        $('#selected-steps').sortable('disable');  
    }
}

async function submitSteps() {
    try {
        collectStepInputs(steps);
        const response = await fetch(`/batch_code/submit_form_code/${task_id}`, {
            method: 'POST', headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ steps})});
            console.log(steps)
        if (response.ok) window.location.href = `/batch_code/submit_form_code/${task_id}`;
    } catch (error) {console.error('Error submitting steps:', error);}
}

