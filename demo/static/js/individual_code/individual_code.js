document.getElementById('WGSBtn').addEventListener('click', function() {initialPreset("WGS", fqType)});


async function initialPreset(presetType, fq_type) {
    try {
        const response = await fetch('/batch_code/get_preset', {
            method: 'POST', headers: {'Content-Type': 'application/json'},
            body: JSON.stringify({presetType: presetType,  mode: fq_type})});
        const data = await response.json();  // 这里返回soft_step列表
        data.forEach(element => {addStep(element.type)});
        initOptionsUI(data)
     	disableAllButtons()
    } catch (error) { console.error('Error initialPreset:', error);}
}


document.getElementById('confirmBtn').addEventListener('click', confirmSteps);
async function confirmSteps() {
    confirmBtn.disabled = true;  // 禁用确认按钮
    for (let i = 0; i < steps.length; i++) {steps[i].id = i;}
    steps[0].selectableOptions=filePaths
    console.log(steps)
    try {
        const response = await fetch('/individual_code/init_div_steps', {
            method: 'POST', headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ steps :steps })});
        const data = await response.json();
        initOptionsUI(data);
        soft_step = data;
    } catch (error) { console.error('Error:', error); }
    
    $('#selected-steps').sortable('disable'); // 禁用左侧按钮卡排序
    document.getElementById('submitBtn').disabled = false; // 允许提交
}


document.getElementById('submitBtn').addEventListener('click', submitSteps);
document.getElementById('unlockBtn').addEventListener('click', function() {initialPreset("WGS");});


async function submitSteps() {
    try {
        soft_step.forEach((step, stepId) => {
            step.params = document.getElementById(`input-${stepId}`).value

            const selectElements = document.querySelectorAll(`#input-step-${stepId} select`);
            selectElements.forEach(select => {
                step.select_option[select.title] = select.value;
            });

            const indexElements = document.querySelectorAll(`#index-step-${stepId} select`);
            indexElements.forEach(select => {
                step.select_index[select.title] = select.value;
            });

        })
        const response = await fetch('/batch_code/store_soft_step', {
            method: 'POST', headers: {'Content-Type': 'application/json'},
            body: JSON.stringify({soft_step: soft_step})
        });
        if (response.ok) {window.location.href = '/batch_code/submit_form_code/' + task_id;}

    }catch (error) {console.error('Error submitting steps:', error);}
}


