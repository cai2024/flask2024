document.getElementById('confirmBtn').addEventListener('click', confirmSteps);
async function update_my_options(stepId) {
    const step = steps.find(s => s.id === stepId);
    try {
        const response = await fetch('/update_my_options', {
            method: 'POST', headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ step_id: step.id, step_type: step.type, select_option: step.select_option, mode: step.mode, temp_options: temp_options })
        });
        return  (await response.json()).temp_options;
    } catch (error) { console.error('Error:', error); }
}


async function index_select(stepId) {
    const step = steps.find(s => s.id === stepId);
    try {
        const response = await fetch('/index_select', {
            method: 'POST', headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ step_type: step.type })
        });
        return (await response.json()).index_options;
    } catch (error) { console.error('Error:', error); }
}


// 初始化
async function confirmSteps() {
    confirmBtn.disabled = true;  // 禁用按钮
    console.log('steps:', steps);
    for (let i = 0; i < steps.length; i++) {
        steps[i].id = i;
    }
    for (const step of steps) {
        await initOptionsUI(step.id);
    }
    console.log('Before disable:', $('#selected-steps'));
    $('#selected-steps').sortable('disable');
    const submitBtn = document.getElementById('submitBtn');
    console.log('Before disable:', $('#selected-steps'));
    submitBtn.disabled = false; 
}


// 初始化当前行的input和options
// /get_options /index_select  返回选项和参数





async function init_type_div(stepId,options) {
    const step = steps.find(s => s.id === stepId);
    try {
        const response = await fetch('/init_type_div', {
            method: 'POST', headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ step_id:stepId ,step_type: step.type, step_params: step.params,mode:step.mode,options:options })
        });
        return (await response.json());
    } catch (error) { console.error('Error:', error); }
}




async function initOptionsUI(stepId) {
    const step = steps.find(step => step.id === stepId);
    const temp_data = await init_type_div(stepId,temp_options)
    const inputDiv = document.createElement('div');
    inputDiv.innerHTML = temp_data.div_string
    document.getElementById('inputArea').appendChild(inputDiv);

    const inputElement = document.querySelector(`#input-step-${stepId} input`);
    step.params = inputElement.value;
    inputElement.addEventListener('input', () => { updateOutputArea(); });

    for (let key in temp_data.rec_option) {
        const selectElement = document.querySelector(`#input-step-${stepId} select#select-${stepId}-${key}`)
        const index = Array.from(selectElement.options).findIndex(option => option.value === temp_data.rec_option[key]);
        if (index !== -1) { selectElement.selectedIndex = index; }
        step.select_option[key] = selectElement.value; 
        temp_options = await update_my_options(stepId);   
    }
    if (step.type.slice(0, 3) === 'map') {
        const index_options =  await index_select(stepId);
        const optionsHtml = index_options.map(option => `<option value="${option}">${option}</option>`).join('');
        const inputDiv = document.createElement('div');
        inputDiv.innerHTML = `
            <div id="index-step-${stepId}" class="form-group">
                <label><strong>step${stepId}: ${step.type} 参考基因组:</strong></label>
                <select class="form-control mt-2">${optionsHtml}</select>
            </div> `;
        document.getElementById('indexArea').appendChild(inputDiv);
        const selectElement = inputDiv.querySelector('select');
        selectElement.onchange = () => { step.select_index=selectElement.value; };
        if ( index_options.length > 0) { step.select_index = selectElement.value; }
    }
}

