$(document).ready(function() {
    const buttons = document.querySelectorAll('.step-button');
    buttons.forEach(button => {
        button.addEventListener('click', function() {
            const buttonType = this.textContent.trim();
            addStep(buttonType, document.querySelector('input[name="fqType"]:checked').value, document.querySelector('input[name="fileType"]:checked').value);
        });
    });

    document.getElementById('WGSBtn').addEventListener('click', function() {initialPreset("WGS",document.querySelector('input[name="fqType"]:checked').value)});
    document.getElementById('confirmBtn').addEventListener('click', confirmSteps);
    document.getElementById('submitBtn').addEventListener('click', submitSteps);
    document.getElementById('unlockBtn').addEventListener('click', function() {initialPreset("WGS");});
});

// 定义了initialPreset， confirmSteps， submitSteps
// // initialPreset读入presetType和mode，addStep获取步骤，initOptionsUI+steps初始化
// // confirmSteps读入steps，initOptionsUI+steps初始化
// // steps获取网页的params,select_option,select_index,上传step+taks_id到session并跳转
// // 审阅时间：20241021


async function initialPreset(presetType, mode) {
    try {
        const response = await fetch('/batch_code/get_preset', {
            method: 'POST', headers: {'Content-Type': 'application/json'},
            body: JSON.stringify({presetType: presetType,  mode: mode})});
        const data = await response.json();
        data.forEach(element => {addStep(element.type, element.mode, element.raw_type)});
        steps = data
        initOptionsUI(data)
    } catch (error) { console.error('Error initialPreset:', error);}
    disableAllButtons()
    confirmBtn.disabled = true;
    $('#selected-steps').sortable('disable');
    document.getElementById('submitBtn').disabled = false;
}

async function confirmSteps() {
    confirmBtn.disabled = true;  // 禁用确认按钮
    for (let i = 0; i < steps.length; i++) {steps[i].id = i;}
    try {
        const response = await fetch('/batch_code/init_div_steps', {
            method: 'POST', headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ steps :steps })});
        const data = await response.json();
        steps = data
        initOptionsUI(data);
    } catch (error) { console.error('Error:', error); }
    $('#selected-steps').sortable('disable'); // 禁用左侧按钮卡排序
    document.getElementById('submitBtn').disabled = false; // 允许提交
}

async function submitSteps() {
    try {
        steps.forEach((step, stepId) => {
            step.params = document.getElementById(`input-${stepId}`).value
            const selectElements = document.querySelectorAll(`#input-step-${stepId} select`);
            selectElements.forEach(select => {step.select_option[select.title] = select.value;});
            const indexElements = document.querySelectorAll(`#index-step-${stepId} select`);
            indexElements.forEach(select => {step.select_index[select.title] = select.value;});
        })
        const response = await fetch('/batch_code/submit_form_code/${task_id}', {
            method: 'POST', headers: {'Content-Type': 'application/json'},
            body: JSON.stringify({steps: steps,task_id: task_id})
        });
        if (response.ok) {window.location.href = '/batch_code/submit_form_code/' + task_id;}
    }catch (error) {console.error('Error submitting steps:', error);}
}

